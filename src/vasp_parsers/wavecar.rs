#![allow(non_upper_case_globals)]

use std::{
    io::{
        Seek,
        SeekFrom,
    },
    fmt,
    fs::File
};

use byteorder::{
    ReadBytesExt,
    LittleEndian,
};
use ndarray::{
    Array1,
    Array2,
    Array3,
    arr1,
    arr2,
    
};
use anyhow::bail;
use ndrustfft::Complex;
//use log::warn;
//use rayon::prelude::*;

use crate::{
    types::{
        Axis,
        Mat33,
        MatX3,
        Result,
    },
    vasp_parsers::binary_io::ReadArray,
};


// Constants
const PI:           f64 = std::f64::consts::PI;
const PIx2:         f64 = PI * 2.0;
const H_PLANCK:     f64 = 6.6260755E-34;
const HBAR:         f64 = H_PLANCK / PIx2;
const RY_TO_EV:     f64 = 13.605693009;
const AU_TO_A:      f64 = 0.529177249;
const HBAR2D2ME:    f64 = RY_TO_EV * AU_TO_A * AU_TO_A;


/// Wavefunction precision type
///
/// Usually VASP stores the band coefficients in two precisions: float32 and float64.
/// From the precision tag in WAVECAR's header we can infer the precision type.
#[derive(PartialEq, Debug, Clone, Copy)]
pub enum WFPrecType {
    Complex32,
    Complex64,
}

/// WAVECAR type enumeration
#[derive(PartialEq, Debug, Clone, Copy)]
pub enum WavecarType {
    /// Most typical WAVECAR type. VASP is executed via `vasp_std`, calculation is done
    /// at all k-points.
    Standard,
    /// VASP is executed via `vasp_gam`, only half of the k-space is utilized because
    /// only Gamma point is sampled, which meas the wavefunction in real space is totally
    /// real, thus only half of the coefficients are needed in k-space.
    GamaHalf(Axis),
    /// VASP is executed via `vasp_ncl`, there should be `LNONCOLLINEAR = T` in OUTCAR.
    /// In this type of WAVECAR, ISPIN = 1, but two spinor components are stored in ONE
    /// spin component side by side.
    NonCollinear,
}

impl fmt::Display for WavecarType {
    fn fmt(&self, f: &mut fmt::Formatter<'_>) -> fmt::Result {
        let description = match self {
            Self::Standard          => "Standard",
            Self::NonCollinear      => "NonCollinear",
            Self::GamaHalf(Axis::X) => "GammaX",
            Self::GamaHalf(Axis::Y) => "GammaY",
            Self::GamaHalf(Axis::Z) => "GammaZ",
        };
        f.write_str(description)
    }
}


pub struct Wavecar {
    file:           File,

    file_len:       u64,
    rec_len:        u64,
    prec_type:      WFPrecType,
    wavecar_type:   WavecarType,

    nspin:          u64,
    nkpoints:       u64,
    nbands:         u64,

    encut:          f64,
    efermi:         f64,

    acell:          Mat33<f64>,
    bcell:          Mat33<f64>,

    volume:         f64,
    ngrid:          [u64; 3],

    nplws:          Vec<u64>,
    kvecs:          Array2<f64>,
    band_eigs:      Array3<f64>,
    band_fweights:  Array3<f64>,
}


impl Wavecar {
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
            ikpoint + (nbands + 1) +
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
            WavecarType::GamaHalf(axis) => Self::_filter_fft_grid(gvecs, axis)
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
            x if x == gvecs_x.len()     => Ok(WavecarType::GamaHalf(Axis::X)),
            x if x == gvecs_z.len()     => Ok(WavecarType::GamaHalf(Axis::Z)),
            _ => bail!("Unknown wavecar type found."),
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

    pub fn read_wavefunction_coeffs(&mut self,
                                    ispin: u64,
                                    ikpoint: u64,
                                    iband: u64) -> Result<Array1<Complex<f64>>> {
        let seek_pos = self.calc_record_location(ispin, ikpoint, iband)?;
        self.file.seek(seek_pos).unwrap();

        let nplw = self.nplws[ikpoint as usize] as usize;
        let dump = match self.prec_type {
            WFPrecType::Complex32 => {
                let mut ret = vec![0f32; nplw * 2];
                self.file.read_f32_into::<LittleEndian>(&mut ret)?;
                ret.into_iter()
                    .map(|x| x as f64)
                    .collect::<Vec<f64>>()
            },
            WFPrecType::Complex64 => {
                let mut ret = vec![0f64; nplw * 2];
                self.file.read_f64_into::<LittleEndian>(&mut ret)?;
                ret
            }
        };

        Ok(
            dump.chunks_exact(2)
                .map(|v| Complex::<f64>::new(v[0], v[1]))
                .collect()
        )
    }
}


#[cfg(test)]
mod tests {
    use super::*;
    use std::cmp::Ordering;
    use ndarray::arr1;

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

}
