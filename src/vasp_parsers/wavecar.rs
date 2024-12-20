#![allow(non_upper_case_globals)]

use std::{
    mem,
    slice,
    ops::Range,
    io::{
        Read,
        Seek,
        SeekFrom,
    },
    fmt::{
        self,
        Write as _,
    },
    fs::File,
    path::Path,
    sync::{
        Arc,
        Mutex,
    },
};

use byteorder::{
    ReadBytesExt,
    LittleEndian,
};
use ndarray::{
    self as nd,
    Array1,
    Array2,
    Array3,
    Array4,
    arr1,
    arr2,
    s,
};
use cauchy::Scalar;
use anyhow::{bail, ensure};
use ndrustfft::{
    FftNum,
    FftHandler,
    R2cFftHandler,
    ndifft,
    ndifft_r2c,
    Complex,
};
use hdf5::File as H5File;

use crate::{
    types::{
        Axis,
        Mat33,
        MatX3,
        Result,
    },
    vasp_parsers::poscar::Poscar,
};


// Constants
#[allow(clippy::approx_constant)]
const PI:           f64 = 3.141_592_653_589_793;
const PIx2:         f64 = PI * 2.0;
//const H_PLANCK:     f64 = 6.6260755E-34;
//const HBAR:         f64 = H_PLANCK / PIx2;
const HBAR:         f64 = 0.6582119281559802;               // eV * fs
const RY_TO_EV:     f64 = 13.605826;
const AU_TO_A:      f64 = 0.529177249;
const AU_TO_DEBYE:  f64 = 2.541746;
const HBAR2D2ME:    f64 = RY_TO_EV * AU_TO_A * AU_TO_A;


#[allow(non_camel_case_types)]
type c64 = Complex<f64>;
#[allow(non_camel_case_types)]
type c32 = Complex<f32>;



// impl Norm for ArrayBase to reduce binary size.
// copied from https://docs.rs/ndarray-linalg/latest/src/ndarray_linalg/norm.rs.html
trait Norm {
    type Output;
    fn norm(&self) -> Self::Output;
}

impl<A, S, D> Norm for nd::ArrayBase<S, D>
where A: Scalar,
      S: nd::Data<Elem=A>,
      D: nd::Dimension,
{
    type Output = A::Real;
    fn norm(&self) -> Self::Output {
        self.iter().map(|x| x.square()).sum::<A::Real>().sqrt()
    }
}



// Wavefunction precision type
//
// Usually VASP stores the band coefficients in two precisions: float32 and float64.
// From the precision tag in WAVECAR's header we can infer the precision type.
#[derive(PartialEq, Debug, Clone, Copy)]
pub enum WFPrecType {
    Complex32,
    Complex64,
}

/// WAVECAR type enumeration
#[derive(PartialEq, Debug, Clone, Copy)]
pub enum WavecarType {
    // Most typical WAVECAR type. VASP is executed via `vasp_std`, calculation is done
    // at all k-points.
    Standard,
    // VASP is executed via `vasp_gam`, only half of the k-space is utilized because
    // only Gamma point is sampled, which meas the wavefunction in real space is totally
    // real, thus only half of the coefficients are needed in k-space.
    GammaHalf(Axis),
    // VASP is executed via `vasp_ncl`, there should be `LNONCOLLINEAR = T` in OUTCAR.
    // In this type of WAVECAR, ISPIN = 1, but two spinor components are stored in ONE
    // spin component side by side.
    NonCollinear,
}

impl fmt::Display for WavecarType {
    fn fmt(&self, f: &mut fmt::Formatter<'_>) -> fmt::Result {
        let description = match self {
            Self::Standard          => "Standard",
            Self::NonCollinear      => "NonCollinear",
            Self::GammaHalf(Axis::X) => "GammaX",
            Self::GammaHalf(Axis::Y) => "GammaY",
            Self::GammaHalf(Axis::Z) => "GammaZ",
        };
        f.write_str(description)
    }
}


#[derive(Clone, Debug)]
pub enum Wavefunction {
    Complex32Array1(Array1<c32>),
    Complex64Array1(Array1<c64>),
    Complex64Array3(Array3<c64>),
    Float64Array3(Array3<f64>),
    Ncl32Array2(Array2<c32>),
    Ncl64Array2(Array2<c64>),
    Ncl64Array4(Array4<c64>),
}


impl Wavefunction {
    pub fn normalize(self) -> Self {
        match self {
            Self::Complex32Array1(wav) => {
                let norm = wav.norm();
                Self::Complex32Array1(wav / norm)
            },
            Self::Complex64Array1(wav) => {
                let norm = wav.norm();
                Self::Complex64Array1(wav / norm)
            },
            Self::Complex64Array3(wav) => {
                let norm = wav.norm();
                Self::Complex64Array3(wav / norm)
            },
            Self::Ncl32Array2(wav) => {
                let norm = wav.norm();
                Self::Ncl32Array2(wav / norm)
            },
            Self::Ncl64Array2(wav) => {
                let norm = wav.norm();
                Self::Ncl64Array2(wav / norm)
            },
            Self::Ncl64Array4(wav) => {
                let norm = wav.norm();
                Self::Ncl64Array4(wav / norm)
            },
            Self::Float64Array3(wav) => {
                let norm = wav.norm();
                Self::Float64Array3(wav / norm)
            },
        }
    }
}


#[derive(Debug)]
pub struct Wavecar {
    file:               Arc<Mutex<File>>,

    pub file_len:       u64,
    pub rec_len:        u64,
    pub prec_type:      WFPrecType,
    pub wavecar_type:   WavecarType,

    pub nspin:          u64,
    pub nkpoints:       u64,
    pub nbands:         u64,

    pub encut:          f64,
    pub efermi:         f64,

    pub acell:          Mat33<f64>,
    pub bcell:          Mat33<f64>,

    pub volume:         f64,
    pub ngrid:          [u64; 3],

    pub nplws:          Vec<u64>,
    pub kvecs:          Array2<f64>,
    pub band_eigs:      Array3<f64>,
    pub band_fweights:  Array3<f64>,
}


impl Wavecar {
    pub fn from_file(path: &(impl AsRef<Path> + ?Sized)) -> Result<Self> {
        let mut file = File::open(path)?;
        let file_len = file.metadata()?.len();

        // Read RECLEN NSPIN and PRECTAG
        file.seek(SeekFrom::Start(0))?;
        let mut dump = [0f64; 3];
        file.read_f64_into::<LittleEndian>(&mut dump)?;
        let rec_len     = dump[0] as u64;
        let nspin       = dump[1] as u64;
        let prec_tag    = dump[2] as u64;

        let prec_type = match prec_tag {
            45200 => WFPrecType::Complex32,
            45210 => WFPrecType::Complex64,
            53300 => bail!("Unsupported WAVECAR: VASP5 with f32."),
            53310 => bail!("Unsupported WAVECAR: VASP5 with f64."),
            _     => bail!("Unknown WAVECAR format."),
        };

        // Read NKPTS NBANDS ENCUT and lattice info
        file.seek(SeekFrom::Start(rec_len))?;
        let mut dump = [0f64; 3 + 9 + 1];
        file.read_f64_into::<LittleEndian>(&mut dump)?;
        let nkpoints = dump[0] as u64;
        let nbands   = dump[1] as u64;
        let encut    = dump[2];

        let mut acell = Mat33::<f64>::default();
        acell.iter_mut().flatten()
            .zip(dump[3 .. 12].iter())
            .for_each(|(ac, du)| {
                *ac = *du;
            });
        let efermi   = dump[12];

        let bcell    = Poscar::mat33_transpose(&Poscar::mat33_inv(&acell).unwrap());
        let volume   = Poscar::mat33_det(&acell);

        if volume <= 1E-5 {
            bail!("WAVECAR corruption: lattice info not correct: ACELL = {:#?} , VOLUME = {}", acell, volume);
        }

        let ngrid = {
            let tmp = acell.iter()
                .map(|[x, y, z]| f64::sqrt(x*x + y*y + z*z))    // lattice vector length
                .map(|len| {
                    ((encut / RY_TO_EV).sqrt() / 
                     (PIx2 / (len / AU_TO_A))).ceil() as u64
                })
                .map(|x| x*2 + 1)
                .collect::<Vec<u64>>();
            [tmp[0], tmp[1], tmp[2]]
        };

        let (nplws, kvecs, band_eigs, band_fweights) =
            Self::_read_band_info(&mut file, nspin, nkpoints, nbands, rec_len)?;

        let wavecar_type = Self::_determine_wavecar_type(&ngrid, &kvecs.row(0).to_owned(), &bcell, encut, nplws[0])?;

        Ok(Self {
            file: Arc::new(Mutex::new(file)),

            file_len,
            rec_len,
            prec_type,
            wavecar_type,

            nspin,
            nkpoints,
            nbands,

            encut,
            efermi,
            acell,
            bcell,

            volume,
            ngrid,

            nplws,
            kvecs,
            band_eigs,
            band_fweights,
        })
    }


    pub fn set_wavecar_type(&mut self, t: WavecarType) -> Result<()> {
        self.check_wavecar_type(t)?;
        self.wavecar_type = t;
        Ok(())
    }


    #[allow(clippy::type_complexity)]
    fn _read_band_info(file:   &mut File,
                       nspin:       u64,
                       nkpoints:    u64,
                       nbands:      u64,
                       rec_len:     u64) -> Result<(Vec<u64>,
                                                    Array2<f64>,
                                                    Array3<f64>,
                                                    Array3<f64>)> {
        let mut nplws           = Vec::<u64>::new();
        let mut kvecs           = Vec::<f64>::new();
        let mut band_eigs       = Vec::<f64>::new();
        let mut band_fweights   = Vec::<f64>::new();

        for ispin in 0 .. nspin {
            for ikpoint in 0 .. nkpoints {
                let rec_idx = Self::_calc_record_index(ispin, ikpoint, 0, nkpoints, nbands);
                let rec_loc = SeekFrom::Start((rec_idx - 1) * rec_len);

                let mut dump = vec![0f64; 4 + 3 * nbands as usize];
                file.seek(rec_loc)?;
                file.read_f64_into::<LittleEndian>(&mut dump)?;

                if 0 == ispin {
                    nplws.push(dump[0] as u64);
                    kvecs.extend_from_slice(&dump[1..4]);
                }

                let dump = dump[4..].to_vec();
                band_eigs.extend(dump.iter().step_by(3));
                band_fweights.extend(dump[2..].iter().step_by(3));
            }
        }

        let kvecs = Array2::from_shape_vec((nkpoints as usize, 3), kvecs)?;
        let band_eigs = Array3::from_shape_vec((nspin as usize, nkpoints as usize, nbands as usize), band_eigs)?;
        let band_fweights = Array3::from_shape_vec((nspin as usize, nkpoints as usize, nbands as usize), band_fweights)?;

        Ok((nplws, kvecs, band_eigs, band_fweights))
    }


    /// Indices starts from 0.
    ///
    /// First 2 means two records for headers
    #[inline]
    fn _calc_record_index(ispin:    u64,
                          ikpoint:  u64,
                          iband:    u64,
                          nkpoints: u64,
                          nbands:   u64) -> u64 {
        2 + ispin * nkpoints * (nbands + 1) +
            ikpoint * (nbands + 1) +
            iband + 1
    }


    #[inline]
    fn _calc_record_location(ispin: u64,
                             ikpoint: u64,
                             iband: u64,
                             nkpoints: u64,
                             nbands: u64,
                             rec_len: u64) -> SeekFrom {
        SeekFrom::Start(
            Self::_calc_record_index(ispin, ikpoint, iband, nkpoints, nbands) * rec_len
        )
    }


    /// Indices starts from 0
    #[inline]
    fn calc_record_location(&self, ispin: u64, ikpoint: u64, iband: u64) -> Result<SeekFrom> {
        self.check_indices(ispin, ikpoint, iband)?;
        Ok(Self::_calc_record_location(ispin, ikpoint, iband, self.nkpoints, self.nbands, self.rec_len))
    }


    /// Indices starts from 0
    #[inline]
    pub fn check_indices(&self, ispin: u64, ikpoint: u64, iband: u64) -> Result<()> {
        self.check_spin_index(ispin)
            .and_then(|_| self.check_kpoint_index(ikpoint))
            .and_then(|_| self.check_band_index(iband))
    }

    #[inline]
    pub fn check_spin_index(&self, ispin: u64) -> Result<()> {
        if ispin >= self.nspin { bail!("Spin index outbound."); }
        Ok(())
    }


    #[inline]
    pub fn check_kpoint_index(&self, ikpoint: u64) -> Result<()> {
        if ikpoint >= self.nkpoints { bail!("K point index outbound."); }
        Ok(())
    }


    #[inline]
    pub fn check_band_index(&self, iband: u64) -> Result<()> {
        if iband >= self.nbands { bail!("Band index outbound."); }
        Ok(())
    }


    fn _generate_fft_freq(ngrid: u64) -> Vec<i64> {
        // ret = [0 ..= ngrid/2] ++ [(1 + ngrid/2 - ngrid) ..= -1];
        // eg: ngrid = 11, ret = vec![0, 1, 2, 3, 4, 5, -5, -4, -3, -2, -1];
        let mut ret = Vec::<i64>::new();
        let ngrid = ngrid as i64;
        ret.extend(
            (0 ..= (ngrid/2)).chain((1+ngrid/2-ngrid) ..= -1)
            );
        ret
    }


    fn _generate_fft_grid_general(ngrid:    &[u64; 3],
                                  kvec:     &Array1<f64>,
                                  bcell:    &Mat33<f64>,
                                  encut:    f64) -> MatX3<i64> {
        let bcell = arr2(bcell);
        let bcell_t = bcell.t();
        let fx = Self::_generate_fft_freq(ngrid[0]);
        let fy = Self::_generate_fft_freq(ngrid[1]);
        let fz = Self::_generate_fft_freq(ngrid[2]);

        let mut ret = vec![[0i64; 3]; 0];

        for z in fz {
            for y in fy.iter().copied() {
                for x in fx.iter().copied() {
                    // G + k
                    let gpk = arr1(&[x as f64 + kvec[0],
                                     y as f64 + kvec[1],
                                     z as f64 + kvec[2]]);
                    let normsqr = {
                        let r = bcell_t.dot(&gpk);
                        r[0] * r[0] + r[1] * r[1] + r[2] * r[2]
                    };
                    if normsqr * PIx2.powi(2) * HBAR2D2ME < encut {
                        ret.push([x, y, z]);
                    }
                }
            }
        }

        ret
    }


    fn _filter_fft_grid(gvecs: MatX3<i64>, axis: Axis) -> MatX3<i64> {
        match axis {
            Axis::X => {
                gvecs.into_iter().filter(|[gx, gy, gz]| {
                    (*gx > 0) ||
                        (*gx == 0 && *gy > 0) ||
                        (*gx == 0 && *gy == 0 && *gz >= 0)
                })
                .collect()
            },

            Axis::Y => {
                unimplemented!("No y-direction gamma half fft scheme was implemented in VASP.");
            },
            
            Axis::Z => {
                gvecs.into_iter().filter(|[gx, gy, gz]| {
                    (*gz > 0) ||
                        (*gz == 0 && *gy > 0) ||
                        (*gz == 0 && *gy == 0 && *gx >= 0)
                })
                .collect()
            }
        }
    }


    fn _generate_fft_grid_specific(ngrid:   &[u64; 3],
                                   kvec:    &Array1<f64>,
                                   bcell:   &Mat33<f64>,
                                   encut:   f64,
                                   wavecar_type: WavecarType) -> MatX3<i64> {
        let gvecs = Self::_generate_fft_grid_general(ngrid, kvec, bcell, encut);
        match wavecar_type {
            WavecarType::Standard | WavecarType::NonCollinear 
                                        => gvecs,
            WavecarType::GammaHalf(axis) => Self::_filter_fft_grid(gvecs, axis)
        }
    }


    /// Returns the kgrid indices corresponds to coefficients in WAVECAR.
    pub fn generate_fft_grid(&self, ikpoint: u64) -> MatX3<i64> {
        Self::_generate_fft_grid_specific(
            &self.ngrid,
            &self.kvecs.row(ikpoint as usize).to_owned(),
            &self.bcell,
            self.encut,
            self.wavecar_type,
        )
    }


    /// Generate kinetic vetor of each FFT grid point, in eV*fs/Angstrom
    pub fn generate_fft_grid_cart(&self, ikpoint: u64) -> MatX3<f64> {
        let kvec  = self.kvecs.row(ikpoint as usize);
        self.generate_fft_grid(ikpoint)
            .into_iter()
            .map(|[gx, gy, gz]| [gx as f64 + kvec[0], gy as f64 + kvec[1], gz as f64 + kvec[2]])
            .map(|[gx, gy, gz]| {   // now in 1/Angstrom
                [
                    gx * self.bcell[0][0] + gy * self.bcell[1][0] + gz * self.bcell[2][0],
                    gx * self.bcell[0][1] + gy * self.bcell[1][1] + gz * self.bcell[2][1],
                    gx * self.bcell[0][2] + gy * self.bcell[1][2] + gz * self.bcell[2][2],
                ]
            })
            .map(|[px, py, pz]| {   // times 2pi and hbar, in eV*fs/Angstrom
                [
                    px * PIx2 * HBAR,
                    py * PIx2 * HBAR,
                    pz * PIx2 * HBAR,
                ]
            })
            .collect::<Vec<_>>()
    }


    fn _determine_wavecar_type(ngrid:   &[u64; 3],
                               kvec:    &Array1<f64>,
                               bcell:   &Mat33<f64>,
                               encut:   f64,
                               nplw:    u64) -> Result<WavecarType> {
        let gvecs   = Self::_generate_fft_grid_general(ngrid, kvec, bcell, encut);
        let gvecs_x = Self::_filter_fft_grid(gvecs.clone(), Axis::X);
        let gvecs_z = Self::_filter_fft_grid(gvecs.clone(), Axis::Z);

        let nplw    = nplw as usize;

        match nplw {
            x if x == gvecs.len()       => Ok(WavecarType::Standard),
            x if x == gvecs.len() * 2   => Ok(WavecarType::NonCollinear),
            x if x == gvecs_x.len()     => Ok(WavecarType::GammaHalf(Axis::X)),
            x if x == gvecs_z.len()     => Ok(WavecarType::GammaHalf(Axis::Z)),
            _ => bail!("Unknown wavecar type found: nplw = {} , ngvecs = {}", nplw, gvecs.len()),
        }
    }


    fn _check_wavecar_type(ngrid: &[u64; 3],
                           kvec: &Array1<f64>,
                           bcell: &Mat33<f64>,
                           encut: f64,
                           nplw: u64,
                           t: WavecarType) -> Result<()> {
        let gvecs = Self::_generate_fft_grid_specific(ngrid, kvec, bcell, encut, t);
        let nplw = nplw as usize;

        if gvecs.len() == nplw || (gvecs.len() * 2 == nplw) {
            return Ok(());
        }

        let suggest_type = Self::_determine_wavecar_type(ngrid, kvec, bcell, encut, nplw as u64)?;
        bail!("Unmatched WAVECAR type: {}, suggested: {}", t, suggest_type);
    }


    pub fn check_wavecar_type(&self, t: WavecarType) -> Result<()> {
        Self::_check_wavecar_type(&self.ngrid,
                                  &self.kvecs.row(0).to_owned(),
                                  &self.bcell,
                                  self.encut,
                                  self.nplws[0],
                                  t)
    }


    fn _read_wavefunction_raw<T: FftNum>(&self,
                                         ispin: u64,
                                         ikpoint: u64,
                                         iband: u64) -> Result<Array1<Complex<T>>> {
        let seek_pos = self.calc_record_location(ispin, ikpoint, iband)?;

        let nplw = self.nplws[ikpoint as usize] as usize;
        let size = nplw * mem::size_of::<Complex<T>>();
        let mut ret = Array1::<Complex<T>>::zeros(nplw);

        unsafe {
            let ptr = ret.as_mut_ptr();
            let ret_slice = slice::from_raw_parts_mut(ptr as *mut u8, size);

            let mut file = self.file.lock().unwrap();
            file.seek(seek_pos)?;
            file.read_exact(ret_slice)?;
        }

        Ok(ret)
    }


    /// Indices starts from 0
    pub fn read_wavefunction(&self,
                             ispin: u64,
                             ikpoint: u64,
                             iband: u64) -> Result<Wavefunction> {
        let nplw = self.nplws[ikpoint as usize] as usize;
        match self.prec_type {
            WFPrecType::Complex32 => {
                let ret = self._read_wavefunction_raw(ispin, ikpoint, iband)?;
                if self.wavecar_type != WavecarType::NonCollinear {     // std & gam wavefunction
                    return Ok(Wavefunction::Complex32Array1(ret));
                }

                ensure!(nplw % 2 == 0, "Odd NPLW for NCL WAVECAR.");   // ncl wavefunction, two spinors
                let nplw = nplw / 2;
                Ok(Wavefunction::Ncl32Array2(ret.into_shape_with_order((2, nplw)).unwrap()))
            },
            WFPrecType::Complex64 => {
                let ret = self._read_wavefunction_raw(ispin, ikpoint, iband)?;
                if self.wavecar_type != WavecarType::NonCollinear {
                    return Ok(Wavefunction::Complex64Array1(ret));
                }

                ensure!(nplw % 2 == 0, "Odd NPLW for NCL WAVECAR.");
                let nplw = nplw / 2;
                Ok(Wavefunction::Ncl64Array2(ret.into_shape_with_order((2, nplw)).unwrap()))
            },
        }
    }


    pub fn show_eigs_fweights(&self) -> String {
        let eigs = self.band_eigs.to_owned() - self.efermi;
        let occs = &self.band_fweights;
        let mut ret = String::new();

        for ispin in 0 .. self.nspin as usize {
            for ikpoint in 0 .. self.nkpoints as usize {
                for iband in 0 .. self.nbands as usize {
                    writeln!(ret, "ISPIN: {:2}  IKPOINT: {:4}  IBAND: {:5} E-Efermi: {:10.3} Occupation: {:5.3}",
                                    ispin + 1, ikpoint + 1, iband + 1,
                                    eigs[[ispin, ikpoint, iband]], occs[[ispin, ikpoint, iband]]).unwrap();
                }
            }
        }

        ret
    }


    // indices start from 0
    pub fn get_wavefunction_realspace(&self, ispin: u64, ikpoint: u64, iband: u64, ngrid: Option<[u64; 3]>) -> Result<Wavefunction> {
        assert!(ispin < self.nspin, "Invalid ispin: {}, nspin = {}", ispin + 1, self.nspin);
        assert!(ikpoint < self.nkpoints, "Invalid ikpoint: {}, nkpoints = {}", ikpoint + 1, self.nkpoints);
        assert!(iband < self.nbands, "Invalid iband: {}, nbands = {}", iband + 1, self.nbands);

        let ngrid = if let Some(g) = ngrid {
            if g[0] < self.ngrid[0] || g[1] < self.ngrid[1] || g[2] < self.ngrid[2] {
                bail!("NGXF NGYF NGZF smaller than NGX NGY NGZ, pleas check.")
            }
            g
        } else {
            [ self.ngrid[0] * 2,
              self.ngrid[1] * 2,
              self.ngrid[2] * 2, ]
        };

        let ngxr = ngrid[0] as i64;
        let ngyr = ngrid[1] as i64;
        let ngzr = ngrid[2] as i64;

        match self.wavecar_type {
            WavecarType::Standard           => self._get_wavefunction_realspace_std(ispin, ikpoint, iband, ngxr, ngyr, ngzr),
            WavecarType::NonCollinear       => self._get_wavefunction_realspace_ncl(ispin, ikpoint, iband, ngxr, ngyr, ngzr),
            WavecarType::GammaHalf(Axis::X)  => self._get_wavefunction_realspace_gamx(ispin, ikpoint, iband, ngxr, ngyr, ngzr),
            WavecarType::GammaHalf(Axis::Z)  => self._get_wavefunction_realspace_gamz(ispin, ikpoint, iband, ngxr, ngyr, ngzr),
            _ => bail!("Unknown or unsupported WAVECAR: {}", self.wavecar_type),
        }
    }


    // All indices counts from 0.
    fn _get_wavefunction_realspace_std(&self, ispin: u64, ikpoint :u64, iband: u64, ngxr: i64, ngyr: i64, ngzr: i64) -> Result<Wavefunction> {
        assert_eq!(self.wavecar_type, WavecarType::Standard);

        let gvecs: MatX3<usize> = self.generate_fft_grid(ikpoint)
            .into_iter()
            .map(|[x, y, z]| {
                let gx = x.rem_euclid(ngxr) as usize;
                let gy = y.rem_euclid(ngyr) as usize;
                let gz = z.rem_euclid(ngzr) as usize;

                [gx, gy, gz]
            })
            .collect();

        let ngxr = ngxr as usize;
        let ngyr = ngyr as usize;
        let ngzr = ngzr as usize;

        let coeffs: Array1<c64> = match self.prec_type {
            WFPrecType::Complex32 => self._read_wavefunction_raw::<f32>(ispin, ikpoint, iband)?
                .mapv(|x| Complex::<f64>{re: x.re as f64, im: x.im as f64}),
            WFPrecType::Complex64 => self._read_wavefunction_raw::<f64>(ispin, ikpoint, iband)?,
        };

        assert_eq!(coeffs.len(), gvecs.len());
        let mut wavk = Array3::<c64>::zeros((ngxr, ngyr, ngzr));
        let mut wavr = Array3::<c64>::zeros((ngxr, ngyr, ngzr));

        gvecs.into_iter().zip(coeffs)
            .for_each(|(idx, v)| wavk[idx] = v);

        let handlers: [FftHandler<f64>; 3] = [
            FftHandler::new(ngxr),
            FftHandler::new(ngyr),
            FftHandler::new(ngzr),
        ];
        ndifft(&wavk, &mut wavr, &handlers[0], 0);
        ndifft(&wavr, &mut wavk, &handlers[1], 1);
        ndifft(&wavk, &mut wavr, &handlers[2], 2);

        Ok(Wavefunction::Complex64Array3(wavr))
    }


    // All indices counts from 0.
    fn _get_wavefunction_realspace_gamx(&self, ispin: u64, ikpoint: u64, iband: u64, ngxr: i64, ngyr: i64, ngzr: i64) -> Result<Wavefunction> {
        assert_eq!(self.wavecar_type, WavecarType::GammaHalf(Axis::X));

        let ngxk = ngxr / 2 + 1;
        let ngyk = ngyr;
        let ngzk = ngzr;

        let gvecs = self.generate_fft_grid(ikpoint)
            .into_iter()
            .map(|[x, y, z]| {
                let gx = x.rem_euclid(ngxk) as usize;
                let gy = y.rem_euclid(ngyk) as usize;
                let gz = z.rem_euclid(ngzk) as usize;

                [gx, gy, gz]
            });

        let ngxr = ngxr as usize;
        let ngyr = ngyr as usize;
        let ngzr = ngzr as usize;

        let ngxk = ngxk as usize;
        let ngyk = ngyk as usize;
        let ngzk = ngzk as usize;

        let coeffs: Array1<c64> = match self.prec_type {
            WFPrecType::Complex32 => self._read_wavefunction_raw::<f32>(ispin, ikpoint, iband)?
                .mapv(|x| Complex::<f64>{re: x.re as f64, im: x.im as f64}),
            WFPrecType::Complex64 => self._read_wavefunction_raw::<f64>(ispin, ikpoint, iband)?,
        };
        let mut wavk = Array3::<c64>::zeros((ngxk, ngyk, ngzk));
        let mut wavr = Array3::<f64>::zeros((ngxr, ngyr, ngzr));

        gvecs.zip(coeffs)
            .for_each(|(idx, v)| wavk[idx] = v);

        {
            let freqs_y = Self::_generate_fft_freq(ngyk as u64);
            let freqs_z = Self::_generate_fft_freq(ngzk as u64);

            for (iy, fy) in freqs_y.iter().cloned().enumerate() {
                for (iz, fz) in freqs_z.iter().cloned().enumerate() {
                    if fy > 0 || (0 == fy && fz >= 0) { continue; }
                    wavk[[0, iy, iz]] = wavk[[0, ngyk-iy-1, ngzk-iz-1]].conj();
                }
            }
        }

        wavk.mapv_inplace(|v| v.unscale(f64::sqrt(2.0)));
        wavk[[0, 0, 0]].scale(f64::sqrt(2.0));

        let mut work = Array3::<c64>::zeros(wavk.dim());
        let handler_x = R2cFftHandler::<f64>::new(ngxr);
        let handler_y =    FftHandler::<f64>::new(ngyr);
        let handler_z =    FftHandler::<f64>::new(ngzr);
        ndifft    (&wavk, &mut work, &handler_y, 1);
        ndifft    (&work, &mut wavk, &handler_z, 2);
        ndifft_r2c(&wavk, &mut wavr, &handler_x, 0);

        Ok(Wavefunction::Float64Array3(wavr))
    }


    // All indices counts from 0.
    fn _get_wavefunction_realspace_gamz(&self, ispin: u64, ikpoint: u64, iband: u64, ngxr: i64, ngyr: i64, ngzr: i64) -> Result<Wavefunction> {
        assert_eq!(self.wavecar_type, WavecarType::GammaHalf(Axis::Z));

        let ngxk = ngxr;
        let ngyk = ngyr;
        let ngzk = ngzr / 2 + 1;

        let gvecs = self.generate_fft_grid(ikpoint)
            .into_iter()
            .map(|[x, y, z]| {
                let gx = x.rem_euclid(ngxk) as usize;
                let gy = y.rem_euclid(ngyk) as usize;
                let gz = z.rem_euclid(ngzk) as usize;

                [gx, gy, gz]
            });

        let ngxr = ngxr as usize;
        let ngyr = ngyr as usize;
        let ngzr = ngzr as usize;

        let ngxk = ngxk as usize;
        let ngyk = ngyk as usize;
        let ngzk = ngzk as usize;

        let coeffs: Array1<c64> = match self.prec_type {
            WFPrecType::Complex32 => self._read_wavefunction_raw::<f32>(ispin, ikpoint, iband)?
                .mapv(|x| Complex::<f64>{re: x.re as f64, im: x.im as f64}),
            WFPrecType::Complex64 => self._read_wavefunction_raw::<f64>(ispin, ikpoint, iband)?
        };

        let mut wavk = Array3::<c64>::zeros((ngxk, ngyk, ngzk));
        let mut wavr = Array3::<f64>::zeros((ngxr, ngyr, ngzr));

        gvecs.zip(coeffs)
            .for_each(|(idx, v)| wavk[idx] = v);

        {
            let freqs_x = Self::_generate_fft_freq(ngxk as u64);
            let freqs_y = Self::_generate_fft_freq(ngyk as u64);

            for (ix, fx) in freqs_x.iter().cloned().enumerate() {
                for (iy, fy) in freqs_y.iter().cloned().enumerate() {
                    if fy > 0 || (0 == fy && fx >= 0) { continue; }
                    wavk[[ix, iy, 0]] = wavk[[ngxk-ix-1, ngyk-iy-1, 0]].conj();
                }
            }
        }

        wavk.mapv_inplace(|v| v.unscale(f64::sqrt(2.0)));
        wavk[[0, 0, 0]].scale(f64::sqrt(2.0));

        let mut work = Array3::<c64>::zeros(wavk.dim());
        let handler_x =    FftHandler::<f64>::new(ngxr);
        let handler_y =    FftHandler::<f64>::new(ngyr);
        let handler_z = R2cFftHandler::<f64>::new(ngzr);
        ndifft    (&wavk, &mut work, &handler_x, 0);
        ndifft    (&work, &mut wavk, &handler_y, 1);
        ndifft_r2c(&wavk, &mut wavr, &handler_z, 2);

        Ok(Wavefunction::Float64Array3(wavr))
    }


    // All indices counts from 0.
    fn _get_wavefunction_realspace_ncl(&self, ispin: u64, ikpoint: u64, iband: u64, ngxr: i64, ngyr: i64, ngzr: i64) -> Result<Wavefunction> {
        assert_eq!(self.wavecar_type, WavecarType::NonCollinear);

        let gvecs: MatX3<usize> = self.generate_fft_grid(ikpoint)
            .into_iter()
            .map(|[x, y, z]| {
                let gx = x.rem_euclid(ngxr) as usize;
                let gy = y.rem_euclid(ngyr) as usize;
                let gz = z.rem_euclid(ngzr) as usize;

                [gx, gy, gz]
            })
            .collect();

        let ngxr = ngxr as usize;
        let ngyr = ngyr as usize;
        let ngzr = ngzr as usize;

        let coeffs: Array1<c64> = match self.prec_type {
            WFPrecType::Complex32 => self._read_wavefunction_raw::<f32>(ispin, ikpoint, iband)?
                .mapv(|x| Complex::<f64>{re: x.re as f64, im: x.im as f64}),
            WFPrecType::Complex64 => self._read_wavefunction_raw::<f64>(ispin, ikpoint, iband)?
        };

        let nplw = self.nplws[ikpoint as usize] as usize / 2;
        let coeffs: Array2<c64> = coeffs.into_shape_with_order((2, nplw)).unwrap();
        let mut wavk = Array4::<c64>::zeros((2, ngxr, ngyr, ngzr));
        let mut wavr = Array4::<c64>::zeros((2, ngxr, ngyr, ngzr));

        for ispinor in 0 .. 2usize {
            let mut wk = wavk.slice_mut(s![ispinor, .., .., ..]);
            let mut wr = wavr.slice_mut(s![ispinor, .., .., ..]);

            gvecs.iter().zip(coeffs.slice(s![ispinor, ..]))
                .for_each(|(idx, v)| wk[*idx] = *v);

            let handlers: [FftHandler<f64>; 3] = [
                FftHandler::new(ngxr),
                FftHandler::new(ngyr),
                FftHandler::new(ngzr),
            ];
            ndifft(&wk, &mut wr, &handlers[0], 0);
            ndifft(&wr, &mut wk, &handlers[1], 1);
            ndifft(&wk, &mut wr, &handlers[2], 2);
        }

        Ok(Wavefunction::Ncl64Array4(wavr))
    }


    /// Read wavefunction coefficients and convert to Array1<Complex64>
    ///
    /// All indices starts from 0.
    pub fn _wav_kspace(&self, ispin: u64, ikpoint: u64, iband: u64, nplw: usize) -> Array2<c64> {
        let wav = self.read_wavefunction(ispin, ikpoint, iband).unwrap();
        match wav {
            Wavefunction::Complex32Array1(wf) => {
                wf.mapv(|x| Complex::<f64>::new(x.re as f64, x.im as f64))
                    .into_shape_with_order((1, nplw))
            },
            Wavefunction::Complex64Array1(wf) => {
                wf.into_shape_with_order((1, nplw))
            },
            Wavefunction::Float64Array3(wf) => {
                wf.mapv(|x| Complex::<f64>::new(x, 0.0))
                    .into_shape_with_order((1, nplw))
            },
            Wavefunction::Ncl32Array2(wf) => {
                wf.mapv(|x| Complex::<f64>::new(x.re as f64, x.im as f64))
                    .into_shape_with_order((2, nplw))
            }
            Wavefunction::Ncl64Array2(wf) => {
                wf.into_shape_with_order((2, nplw))
            }
            _ => panic!("Unexpected wavefunction type passed to '_wav_kspace'.")
        }
        .unwrap()
    }


    /// Calculate transition dipole moment using pseudo wavefunction.
    ///
    /// All indices starts from 0.
    pub fn transition_dipole_moment(&self, ispin: u64, ikpoint: u64, iniband: u64, finband: u64) -> [c64; 3] {
        self.check_indices(ispin, ikpoint, iniband).unwrap();
        self.check_indices(ispin, ikpoint, finband).unwrap();
        assert!(iniband != finband);

        #[allow(non_snake_case)]
        let dE = self.band_eigs[(ispin as usize, ikpoint as usize, finband as usize)] - 
                 self.band_eigs[(ispin as usize, ikpoint as usize, iniband as usize)];

        let kvec = self.kvecs.row(ikpoint as usize);
        let gvecs = self.generate_fft_grid(ikpoint)
            .into_iter()
            .map(|[gx, gy, gz]| [gx as f64 + kvec[0], gy as f64 + kvec[1], gz as f64 + kvec[2]])
            .collect::<Vec<_>>();
        let gvecs = arr2(&gvecs)
            .dot(&(arr2(&self.bcell) * PIx2))
            .mapv(|v| Complex::<f64>::new(v, 0.0));

        let nplw = gvecs.shape()[0];
        let nspinor = if WavecarType::NonCollinear == self.wavecar_type {
            2
        } else {
            1
        };
        assert_eq!(nplw * nspinor, self.nplws[ikpoint as usize] as usize);

        let phi_i = self._wav_kspace(ispin, ikpoint, iniband, nplw);
        let phi_j = self._wav_kspace(ispin, ikpoint, finband, nplw);

        let olap = match self.wavecar_type {
            WavecarType::Standard | WavecarType::NonCollinear => {
                phi_j.mapv_into(|v| v.conj()) * phi_i  // phi_j.conj() .* phi_i
            },
            _ => {
                (phi_j.mapv(|v| v.conj()) * &phi_i -
                 phi_i.mapv(|v| v.conj()) * &phi_j) / 2.0 // phi_j.conj() .* phi_i .- phi_i.conj() .* phi_j
            }
        };

        
        let tdm = (
                olap.dot(&gvecs)    // <phi_j | k | phi_i>
                .sum_axis(ndarray::Axis(0))
                * Complex::<f64>::i()
            )
            .mapv_into(|v| v.scale(AU_TO_A * AU_TO_DEBYE * 2.0 * RY_TO_EV / dE));

        [tdm[0], tdm[1], tdm[2]]
    }


    /// Performs <psi | sigma_z | psi> for given ncl wavefunction.
    pub fn get_sigmaz(psi: &Array2<c64>) -> f64 {
        psi.slice(s![0, ..]).norm() - psi.slice(s![1, ..]).norm()
    }


    /// Performs <psi | sigma_z | psi> for ncl wavefunction.
    ///
    /// This method can is dedicated for the ncl system, thus ispin is bounded to be 0
    pub fn get_band_sigmaz(&self, ikpoint: u64, iband: u64) -> Result<f64> {
        ensure!(self.wavecar_type == WavecarType::NonCollinear);
        let nplw = self.nplws[ikpoint as usize] / 2;
        let wav = self._wav_kspace(0, ikpoint, iband, nplw as usize);

        Ok(Self::get_sigmaz(&wav))
    }


    /// Performs <psi_j | sigma_z | psi_i> for given ncl wavefunction pair. psi_i and psi_j
    /// must have same sizes.
    pub fn get_sigmaz_ji(psi_i: &Array2<c64>, psi_j: &Array2<c64>) -> c64 {
        (psi_j.slice(s![0, ..]).mapv(|v| v.conj()) * psi_i.slice(s![0, ..])).sum()
            -
        (psi_j.slice(s![1, ..]).mapv(|v| v.conj()) * psi_i.slice(s![1, ..])).sum()
    }


    /// Performs <psi_j | sigma_z | psi_i> for ncl wavefunctions. Where psi_i and psi_j
    /// must be in same k-point.
    ///
    /// This method can is dedicated for the ncl system, thus ispin is bounded to be 0
    pub fn get_band_sigmaz_ji(&self, ikpoint: u64, iband: u64, jband: u64) -> Result<c64> {
        ensure!(self.wavecar_type == WavecarType::NonCollinear);
        let nplw = self.nplws[ikpoint as usize] / 2;
        let psi_i = self._wav_kspace(0, ikpoint, iband, nplw as usize);
        let psi_j = self._wav_kspace(0, ikpoint, jband, nplw as usize);

        Ok(
            (psi_j.slice(s![0, ..]).mapv(|v| v.conj()) * psi_i.slice(s![0, ..])).sum()
            -
            (psi_j.slice(s![1, ..]).mapv(|v| v.conj()) * psi_i.slice(s![1, ..])).sum()
        )
    }



    /// Performs an unitary transform for a pair of degenerated state.
    ///
    /// Returns the transformed wavefunction pair where the `sigma_z` of first one purely
    /// points to `+z` and `sigma_z` of the other one points to `-z`.
    ///
    /// NOTE: Please make sure the input wavefunctions are degenerated.
    pub fn pair_unitary_transform(mut psi_i: Array2<c64>, mut psi_j: Array2<c64>, threshold: f64) -> Result<[Array2<c64>; 2]> {
        ensure!(threshold > 0.0);

        let mut z11 = Self::get_sigmaz(&psi_i);
        let mut z22 = Self::get_sigmaz(&psi_j);
        let mut z21 = Self::get_sigmaz_ji(&psi_i, &psi_j);

        if f64::abs(z11 + z22) > threshold {
            bail!("Degeneracy check failed: |z11:{z11:10.5} + z22:{z22:10.5}| = {:10.5} > threshold:{threshold:10.5}", f64::abs(z11+z22));
        }

        let reverse = z11 < 0.0;
        if reverse {
            mem::swap(&mut psi_i, &mut psi_j);
            mem::swap(&mut z11, &mut z22);
            z21.im = -z21.im;     // z21 = z21.conj()
        }

        let expitheta = z21 / z21.norm();
        let phi = f64::atan( z21.norm() / z11 ) / 2.0;

        let a = phi.cos();
        let b = expitheta * phi.sin();

        // U = [[ a , b ],
        //      [-b*, a*]]
        // [psi_i'; psi_j'] = U * [psi_i; psi_j]
        let ret_psi_i = psi_i.clone() * a   + psi_j.clone() * b;
        let ret_psi_j = psi_i * (-b.conj()) + psi_j * a;

        Ok([ret_psi_i, ret_psi_j])
    }


    /// Performs an unitary transform for a pair of degenerated state on same k-point.
    ///
    /// Returns the transformed wavefunction pair where the `sigma_z` of first one purely
    /// points to `+z` and `sigma_z` of the other one points to `-z`.
    ///
    /// NOTE: Please make sure the selected bands are degenerated and have nearly opposite spin
    /// polarization.
    pub fn band_pair_unitary_transform(&self, ikpoint: u64, iband: u64, jband: u64, threshold: f64) -> Result<[Array2<c64>; 2]> {
        ensure!(self.wavecar_type == WavecarType::NonCollinear);
        let nplw = self.nplws[ikpoint as usize] / 2;
        let psi_i = self._wav_kspace(0, ikpoint, iband, nplw as usize);
        let psi_j = self._wav_kspace(0, ikpoint, jband, nplw as usize);
        Self::pair_unitary_transform(psi_i, psi_j, threshold)
    }


    /// Save slices of wavefunction into HDF5 file.
    ///
    /// All indices count from 0.
    pub fn save_slices_to_h5<P: AsRef<Path>>(
        &self, fname: P,
        kslice: Range<usize>,
        bslice: Range<usize>) -> Result<()> {

        let nkpt = kslice.len();
        let nbnd = bslice.len();
        let nspn = self.nspin as usize;
        ensure!(nkpt >= 1 && nkpt < self.nkpoints as usize);
        ensure!(nbnd >= 1 && nbnd < self.nbands as usize);

        let wavtype = match self.wavecar_type {
            WavecarType::Standard => b'S',
            WavecarType::NonCollinear => b'N',
            WavecarType::GammaHalf(Axis::X) => b'X',
            WavecarType::GammaHalf(Axis::Z) => b'Z',
            _ => bail!("Unexpected wavecar type."),
        };

        let prec: usize = match self.prec_type {
            WFPrecType::Complex32 => 32,
            WFPrecType::Complex64 => 64,
        };

        // Range are not Copy, Fxxk off.
        let kid_list = kslice.clone().map(|x| x + 1).collect::<Vec<usize>>();
        let bid_list = bslice.clone().map(|x| x + 1).collect::<Vec<usize>>();

        let nplw_list = self.nplws[kslice.clone()].to_owned();
        let klist = self.kvecs.slice(s![kslice.clone(), ..]).to_owned();
        let eigs = self.band_eigs.slice(s![.., kslice.clone(), bslice.clone()]).to_owned();
        let whts = self.band_eigs.slice(s![.., kslice.clone(), bslice.clone()]).to_owned();

        let f = H5File::open(fname)?;
        f.new_dataset::<usize>().create("prectype")?.write_scalar(&prec)?;
        f.new_dataset::<u8>().create("wavtype")?.write_scalar(&wavtype)?;
        f.new_dataset::<usize>().create("nspin")?.write_scalar(&nspn)?;
        f.new_dataset::<usize>().create("nkpoint")?.write_scalar(&nkpt)?;
        f.new_dataset::<usize>().create("nband")?.write_scalar(&nbnd)?;

        f.new_dataset::<f64>().create("encut")?.write_scalar(&self.encut)?;
        f.new_dataset::<f64>().create("efermi")?.write_scalar(&self.efermi)?;
        f.new_dataset::<[[f64;3];3]>().create("acell")?.write_scalar(&self.acell)?;
        f.new_dataset::<[[f64;3];3]>().create("bcell")?.write_scalar(&self.bcell)?;
        f.new_dataset::<f64>().create("volume")?.write_scalar(&self.volume)?;
        f.new_dataset::<[u64;3]>().create("ngrid")?.write_scalar(&self.ngrid)?;

        f.new_dataset_builder().with_data(&kid_list).create("kpoint_id_list")?;
        f.new_dataset_builder().with_data(&bid_list).create("band_id_list")?;
        f.new_dataset_builder().with_data(&klist).create("kpoint_list")?;
        f.new_dataset_builder().with_data(&eigs).create("band_eigvals")?;
        f.new_dataset_builder().with_data(&whts).create("fermi_weights")?;
        f.new_dataset_builder().with_data(&nplw_list).create("nplws_list")?;


        todo!()
    }
}


impl fmt::Display for Wavecar {
    fn fmt(&self, f: &mut fmt::Formatter<'_>) -> fmt::Result {
        writeln!(f, "{:<10} = {:>20}",      "LENGTH",   self.file_len)?;
        writeln!(f, "{:<10} = {:>20}",      "RECLEN",   self.rec_len)?;
        #[allow(clippy::to_string_in_format_args)]
        writeln!(f, "{:<10} = {:>20}",      "TYPE",     self.wavecar_type.to_string())?;
        writeln!(f, "{:<10} = {:>20}",      "NSPIN",    self.nspin)?;
        writeln!(f, "{:<10} = {:>20}",      "NKPTS",    self.nkpoints)?;
        writeln!(f, "{:<10} = {:>20}",      "NBANDS",   self.nbands)?;
        writeln!(f, "{:<10} = {:>20.3}",    "ENCUT",    self.encut)?;
        writeln!(f, "{:<10} = {:>20.3}",    "EFERMI",   self.efermi)?;
        writeln!(f, "{:<10} = {:>20.3}",    "VOLUME",   self.volume)?;
        writeln!(f, "{:<10} = {:>4?}",      "NGRID",    self.ngrid)?;
        writeln!(f, "{:<10} = {:>5?}",      "NPLWS",    self.nplws)?;
        writeln!(f, "{:<10} = {:>8.3?}",    "ACELL",    self.acell)?;
        writeln!(f, "{:<10} = {:>8.3?}",    "BCELL",    self.bcell)?;
        writeln!(f)?;

        if f.alternate() {
            f.write_str(&self.show_eigs_fweights())?;
        }

        Ok(())
    }
}


#[cfg(test)]
mod tests {
    use super::*;
    use std::cmp::Ordering;
    use ndarray::arr1;

    use crate::vasp_parsers::chg;
    use crate::vasp_parsers::poscar;

    fn generate_fft_grid_ref(ngrid: u64) -> Vec<i64> {
        let ngrid = ngrid as i64;

        (0 .. ngrid).map(|x| match x.cmp(&(ngrid/2 + 1)) {
            Ordering::Less => x,
            _ => x - ngrid,
        })
        .collect() }

    #[test]
    fn test_generate_fft_freq() {
        for i in 2 ..= 50 {
            assert_eq!(Wavecar::_generate_fft_freq(i as u64), generate_fft_grid_ref(i as u64));
        }
    }

    #[test]
    fn test_generate_fft_grid_general() {
        let kvec = arr1(&[1.0/3.0, 1.0/3.0, 0.0]);
        let ngrid = &[11u64, 11, 105];
        let bcell = &[[0.313971743, 0.181271670, 0.000000000],
                      [0.000000000, 0.362543340, 0.000000000],
                      [0.000000000, 0.000000000, 0.028571429]];
        let encut = 323.36125000000004; // little difference from it in OUTCAR: 323.4

        let res = Wavecar::_generate_fft_grid_general(ngrid, &kvec, bcell, encut);
        assert_eq!(res.len(), 3981);
    }

    #[test]
    fn test_determine_wavecar_type() {
        let kvec  = arr1(&[1.0/3.0, 1.0/3.0, 0.0]);
        let ngrid = &[11u64, 11, 105];
        let bcell = &[[0.313971743, 0.181271670, 0.000000000],
                      [0.000000000, 0.362543340, 0.000000000],
                      [0.000000000, 0.000000000, 0.028571429]];
        let encut = 323.36125000000004; // little difference from it in OUTCAR: 323.4

        let wavecar_type = Wavecar::_determine_wavecar_type(ngrid, &kvec, bcell, encut, 3981).unwrap();
        assert_eq!(WavecarType::Standard, wavecar_type);

        let wavecar_type = Wavecar::_determine_wavecar_type(ngrid, &kvec, bcell, encut, 3981 * 2).unwrap();
        assert_eq!(WavecarType::NonCollinear, wavecar_type);
    }

    #[test]
    fn test_calc_record_index() {
        assert_eq!(Wavecar::_calc_record_index(0, 0, 0, 10, 10), 3);
    }

    #[test]
    #[ignore]
    fn test_print_eigs_and_wavefunction() {
        let wav = Wavecar::from_file("WAVECAR").unwrap();
        println!("{}", &wav);
        println!("{:#}", &wav);

        println!("wfc111:\n{:?}", wav.read_wavefunction(0, 0, 0).unwrap());
        println!("wfc222:\n{:?}", wav.read_wavefunction(1, 1, 1).unwrap());
    }

    #[test]
    #[ignore]
    fn test_wavefunction_realspace_std() {
        let wav = Wavecar::from_file("WAVECAR").unwrap();
        let ngxr = wav.ngrid[0] as i64 * 2;
        let ngyr = wav.ngrid[0] as i64 * 2;
        let ngzr = wav.ngrid[0] as i64 * 2;

        for i in 0 .. 1 {
            let wavr = wav._get_wavefunction_realspace_std(0, i, 0, ngxr, ngyr, ngzr).unwrap(); // 1 (i+1) 1
            let mut wavr = match wavr {
                Wavefunction::Complex64Array3(dump) => dump,
                _ => panic!(),
            };

            let normfact = (wavr.len() as f64).sqrt();
            wavr.mapv_inplace(|v| v.scale(normfact));

            println!("{:.10E}\n{:.10E}\n{:.10E}\n", wavr[[0, 0, 0]], wavr[[0, 0, 1]], wavr[[0, 0, 2]]);

            let shape = wavr.shape();

            let chgd = wavr.map(|v| v.re as f64);
            let ngrid = [shape[0], shape[1], shape[2]];

            let pos = poscar::Poscar::from_file("POSCAR").unwrap();
            let chg = chg::ChargeDensity {
                chgtype: chg::ChargeType::Chgcar,
                pos,
                ngrid,
                chg: vec![chgd],
                aug: vec![],
            };
            chg.to_file(&format!("{}.vasp", i)).unwrap();
        }

    }

    #[test]
    #[ignore]
    fn test_wavefunction_realspace_ncl() {
        let wav = Wavecar::from_file("WAVECAR").unwrap();
        let ngxr = wav.ngrid[0] as i64 * 2;
        let ngyr = wav.ngrid[0] as i64 * 2;
        let ngzr = wav.ngrid[0] as i64 * 2;

        for i in 0 .. 64 {
            let wavr = wav._get_wavefunction_realspace_ncl(0, 0, i, ngxr, ngyr, ngzr).unwrap(); // 1 (i+1) 1
            let mut wavr = match wavr {
                Wavefunction::Ncl64Array4(dump) => dump,
                _ => panic!(),
            };

            let normfact = (wavr.len() as f64).sqrt();
            wavr.mapv_inplace(|v| v.scale(normfact));

            println!("{:.10E}\n{:.10E}\n{:.10E}\n", wavr[[0, 0, 0, 0]], wavr[[0, 0, 0, 1]], wavr[[0, 0, 0, 2]]);

            let shape = wavr.shape();

            let chgd = wavr.slice(s![0, .., .., ..]).map(|v| v.re as f64);
            let ngrid = [shape[1], shape[2], shape[3]];

            let pos = poscar::Poscar::from_file("POSCAR").unwrap();
            let chg = chg::ChargeDensity {
                chgtype: chg::ChargeType::Chgcar,
                pos,
                ngrid,
                chg: vec![chgd],
                aug: vec![],
            };
            chg.to_file(&format!("{}.vasp", i)).unwrap();
        }

    }

    #[test]
    #[ignore]
    fn test_wavefunction_realspace_gamx() {
        let wav = Wavecar::from_file("WAVECAR").unwrap();
        let ngxr = wav.ngrid[0] as i64 * 2;
        let ngyr = wav.ngrid[0] as i64 * 2;
        let ngzr = wav.ngrid[0] as i64 * 2;

        for i in 0 .. 1 {
            let wavr = wav._get_wavefunction_realspace_gamx(0, i, 0, ngxr, ngyr, ngzr).unwrap(); // 1 (i+1) 1
            let mut wavr = match wavr {
                Wavefunction::Float64Array3(dump) => dump,
                _ => panic!(),
            };

            let normfact = (wavr.len() as f64).sqrt();
            wavr.mapv_inplace(|v| v / normfact);

            println!("{:.10E}\n{:.10E}\n{:.10E}\n", wavr[[0, 0, 0]], wavr[[0, 0, 1]], wavr[[0, 0, 2]]);

            let shape = wavr.shape();

            let chgd = wavr.map(|v| *v as f64);
            let ngrid = [shape[0], shape[1], shape[2]];

            let pos = poscar::Poscar::from_file("POSCAR").unwrap();
            let chg = chg::ChargeDensity {
                chgtype: chg::ChargeType::Chgcar,
                pos,
                ngrid,
                chg: vec![chgd],
                aug: vec![],
            };
            chg.to_file(&format!("{}.vasp", i)).unwrap();
        }

    }

    #[test]
    #[ignore]
    fn test_wavefunction_realspace_gamz() {
        let mut wav = Wavecar::from_file("WAVECAR").unwrap();
        wav.set_wavecar_type(WavecarType::GammaHalf(Axis::Z)).unwrap();
        let ngxr = wav.ngrid[0] as i64 * 2;
        let ngyr = wav.ngrid[0] as i64 * 2;
        let ngzr = wav.ngrid[0] as i64 * 2;

        for i in 0 .. 1 {
            let wavr = wav._get_wavefunction_realspace_gamz(0, i, 0, ngxr, ngyr, ngzr).unwrap(); // 1 (i+1) 1
            let mut wavr = match wavr {
                Wavefunction::Float64Array3(dump) => dump,
                _ => panic!(),
            };

            let normfact = (wavr.len() as f64).sqrt();
            wavr.mapv_inplace(|v| v / normfact);

            println!("{:.10E}\n{:.10E}\n{:.10E}\n", wavr[[0, 0, 0]], wavr[[0, 0, 1]], wavr[[0, 0, 2]]);

            let shape = wavr.shape();

            let chgd = wavr.mapv(|v| v as f64);
            let ngrid = [shape[0], shape[1], shape[2]];

            let pos = poscar::Poscar::from_file("POSCAR").unwrap();
            let chg = chg::ChargeDensity {
                chgtype: chg::ChargeType::Chgcar,
                pos,
                ngrid,
                chg: vec![chgd],
                aug: vec![],
            };

            chg.to_file(&format!("{}.vasp", i)).unwrap();
        }

    }
}
