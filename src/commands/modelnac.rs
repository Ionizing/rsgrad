use std::path::PathBuf;
use clap::Args;
use anyhow::{
    Context,
    bail,
    ensure,
};
use log::{
    info,
    warn,
};
use ndarray as na;
use hdf5::File as H5File;
use ndrustfft::Complex;

use crate::{
    types::{
        Result,
        Axis,
    },
    OptProcess,
    vasp_parsers::{
        procar::Procar,
        wavecar::{
            Wavecar,
            WavecarType,
            Norm,
        },
        soc::calc_hmm,
    }
};


#[allow(non_camel_case_types)]
type c64 = Complex<f64>;


#[derive(Debug, Args)]
/// Calculate model non-adiabatic coupling (NAC) for NAMD-LMI (a subset of Hefei-NAMD).
///
/// This NAC contains eigenvalues and transition dipole moment (TDM) calculated from selected WAVECAR.
/// The phonon contribution is wipped out, which means the eigenvalues and TDM stay staic over the
/// time, and the NAC (<i|d/dt|j>) vanishes.
///
/// Detailed fields of the produced file:{n}
/// - ikpoint: K point index, counts from 1;{n}
/// - nspin: number of spin channels;{n}
/// - spin_diabatics: use spin diabatics representation or not;{n}
/// - nbands: Number of total bands in WAVECAR;{n}
/// - brange: selected band indices, <start> ..= <end>, index counts from 1;{n}
/// - nbrange: number of selected bands, nbrange = end - start + 1;{n}
/// - efermi: Fermi's level;{n}
/// - normalization: normalize pseudo wavefunction in WAVECAR of not;{n}
/// - potim: ionic time step;{n}
/// - temperature: 1E-6 Kelvin as default;{n}
/// - eigs: band eigenvalues;{n}
/// - pij_r/pij_i: real and imaginary part of <i|p|j>;{n}
/// - proj: Projection on each orbitals of selected bands, cropped from PROCAR.
///
/// Some fields not listed here are not meaningful but essential for the NAMD-LMI.
pub struct ModelNac {
    #[arg(short='w', long, default_value = "./WAVECAR")]
    /// WAVECAR file name.
    wavecar: PathBuf,

    #[arg(long, value_parser = ["x", "z"])]
    /// Gamma Half direction of WAVECAR. You need to set this to 'x' or 'z' when
    /// processing WAVECAR produced by `vasp_gam`.
    gamma_half: Option<String>,

    #[arg(short='k', long, default_value_t = 1)]
    /// One selected K-Point index, count starts from 1.
    ///
    /// Example: --ikpoint 2
    ikpoint: usize,

    #[arg(long, default_value_t = false)]
    /// Use spin diabatics representation or not.
    spin_diabatics: bool,

    #[arg(long, default_value = "./")]
    /// Path that contains NormalCAR and SocCar
    spin_diabatics_path: PathBuf,

    #[arg(long, num_args(2))]
    /// Selected band range, starts from 1.
    ///
    /// Example: --brange 24 42
    brange: Vec<usize>,

    #[arg(long, default_value_t = false)]
    /// Use normalized wavefunctions in WAVECAR to calculate NAC and TDM.
    normalization: bool,

    #[arg(long, default_value_t = 1.0)]
    /// Ionic time step in femtosecond (fs).
    potim: f64,

    #[arg(short='p', long, default_value = "./PROCAR")]
    /// PROCAR file name for band projection parsing.
    procar: PathBuf,

    #[arg(short='o', long, default_value = "./NAC-0K.h5")]
    /// Output file name.
    h5out: PathBuf,
}


impl OptProcess for ModelNac {
    fn process(&self) -> Result<()> {

        info!("Reading WAVECAR: {:?}", &self.wavecar);
        let mut wav = Wavecar::from_file(&self.wavecar)?;
        if let Some(gammahalf) = self.gamma_half.as_ref() {
            if wav.wavecar_type == WavecarType::Standard ||
               wav.wavecar_type == WavecarType::NonCollinear {
                    bail!("Current WAVECAR is not gamma-halved, rsgrad can determine the WAVECAR type directly, \
please remove the argument `gamma_half`.")
            }

            let gammahalf = match gammahalf.as_ref() {
                "x" => WavecarType::GammaHalf(Axis::X),
                "z" => WavecarType::GammaHalf(Axis::Z),
                _ => panic!("Unreachable branch"),
            };
            
            wav.set_wavecar_type(gammahalf)?;
        } else if wav.wavecar_type != WavecarType::Standard &&
            wav.wavecar_type != WavecarType::NonCollinear {
                warn!("Current WAVECAR is gamma-halved, sometimes the gamma-x and gamma-z verions have same plane wave numbers.
I suggest providing `gamma_half` argument to avoid confusion.");
        }
        // cancel mutability
        let wav = wav;

        info!("Reading PROCAR: {:?}", &self.procar);
        let procar = Procar::from_file(&self.procar)?;


        let mut brange = self.brange.clone();
        brange.sort();
        brange.dedup();
        ensure!(brange.len() == 2, "You must input two unique band index.");
        let brange = [brange[0], brange[1]];

        
        let nsw     = 9;
        let nspin   = wav.nspin as usize;
        let spin_diabatics = self.spin_diabatics;
        let nkpoints = wav.nkpoints as usize;
        let ikpoint = self.ikpoint;
        let nbands  = wav.nbands as usize;
        let nbrange = brange[1] - brange[0] + 1;

        let efermi = wav.efermi;
        let olaps = na::Array4::<f64>::zeros((nsw, nspin, nbrange, nbrange));
        let eigs  = wav.band_eigs.slice(na::s![na::NewAxis, .., ikpoint-1, brange[0]-1 .. brange[1]]).to_owned();

        if spin_diabatics && nspin != 2 {
            bail!("Spin diabatics representation requires ISPIN = 2.");
        }

        // <i|p|j>, transition dipole moment
        let mut pijs  = na::Array5::<c64>::zeros((nsw, nspin, 3, nbrange, nbrange));

        let lncl    = match wav.wavecar_type {
            WavecarType::NonCollinear => true,
            _ => false,
        };
        let nspinor = if lncl { 2usize } else { 1 };
        let nplw = wav.nplws[ikpoint - 1] as usize;
        let mut phi = na::Array2::<c64>::zeros((nbrange, nplw));
        let gvecs = na::arr2(&wav.generate_fft_grid_cart(ikpoint as u64 - 1))
            .rows()
            .into_iter()
            .map(|g| [
                c64::new(g[0], 0.0),
                c64::new(g[1], 0.0),
                c64::new(g[2], 0.0),
            ])
            .cycle()
            .take(nplw)
            .flatten()
            .collect::<na::Array1<c64>>()
            .into_shape_with_order((nplw, 3))
            .unwrap();
        for ispin in 0 .. nspin {
            for (ii, iband) in (brange[0] - 1 .. brange[1]).enumerate() {
                phi.slice_mut(na::s![ii, ..]).assign(&{
                    let mut ket = wav._wav_kspace(ispin as u64, ikpoint as u64 - 1, iband as u64, nplw / nspinor)
                        .into_shape_with_order((nplw,))
                        .with_context(|| "Wavefunction reshape failed.")?;
                    if self.normalization {
                        let norm_inv = 1.0 / ket.norm();
                        ket.mapv_inplace(|v| v.scale(norm_inv));
                    }
                    ket
                });
            }

            for idirect in 0 .. 3 {
                let phi_x_gvecs: na::Array2<_> = phi.clone() * gvecs.slice(na::s![na::NewAxis, .., idirect]);

                // <i | p | j>, in eV*fs/Angstrom
                let pij_tmp = match wav.wavecar_type {
                    WavecarType::GammaHalf(_) => phi.mapv(|v| v.conj()).dot(&phi_x_gvecs.t())
                                               - phi_x_gvecs.mapv(|v| v.conj()).dot(&phi.t()),
                    _ => phi.mapv(|v| v.conj()).dot(&phi_x_gvecs.t()),
                };
                pijs.slice_mut(na::s![.., ispin, idirect, .., ..]).assign(&pij_tmp.slice(na::s![na::NewAxis, .., ..]));
            }
        }

        let proj = procar.pdos.projected.slice(na::s![.., ikpoint-1, brange[0]-1 .. brange[1], .., ..]).to_owned();
        let proj = na::stack(na::Axis(0), &vec![proj.view(); nsw])?;

        // Calculate SOC matrix for spin diabatics representation
        let soc = if self.spin_diabatics {
            let hmm = calc_hmm(&self.spin_diabatics_path, nbands, nkpoints, ikpoint)?;
            Some(hmm.slice(na::s![.., brange[0]-1 .. brange[1], brange[0]-1 .. brange[1]]).to_owned())
        } else {
            None
        };

        info!("Saving to {:?}", &self.h5out);
        
        let f = H5File::create(&self.h5out)?;

        f.new_dataset::<usize>().create("ikpoint")?.write_scalar(&self.ikpoint)?;
        f.new_dataset::<usize>().create("nspin")?.write_scalar(&nspin)?;
        f.new_dataset::<bool>().create("spin_diabatics")?.write_scalar(&spin_diabatics)?;
        f.new_dataset::<bool>().create("lncl")?.write_scalar(&lncl)?;
        f.new_dataset::<usize>().create("nbands")?.write_scalar(&nbands)?;
        f.new_dataset::<usize>().create("ndigit")?.write_scalar(&4)?;
        f.new_dataset::<[usize;2]>().create("brange")?.write_scalar(&brange)?;
        f.new_dataset::<usize>().create("nbrange")?.write_scalar(&nbrange)?;
        f.new_dataset::<usize>().create("nsw")?.write_scalar(&(nsw+1))?;
        f.new_dataset::<f64>().create("efermi")?.write_scalar(&efermi)?;
        f.new_dataset::<f64>().create("potim")?.write_scalar(&self.potim)?;
        f.new_dataset::<f64>().create("temperature")?.write_scalar(&1E-6)?;
        f.new_dataset::<bool>().create("normalization")?.write_scalar(&self.normalization)?;
        f.new_dataset::<bool>().create("phasecorrection")?.write_scalar(&true)?;

        f.new_dataset_builder().with_data(&olaps).create("olaps_r")?;
        f.new_dataset_builder().with_data(&olaps).create("olaps_i")?;

        f.new_dataset_builder().with_data(&eigs).create("eigs")?;

        f.new_dataset_builder().with_data(&pijs.mapv(|v| v.re)).create("pij_r")?;
        f.new_dataset_builder().with_data(&pijs.mapv(|v| v.im)).create("pij_i")?;

        if let Some(soc) = soc {
            f.new_dataset_builder().with_data(&soc.mapv(|v| v.re)).create("soc_r")?;
            f.new_dataset_builder().with_data(&soc.mapv(|v| v.im)).create("soc_i")?;
        }
        
        f.new_dataset_builder().with_data(&proj).create("proj")?;

        Ok(())
    }
}
