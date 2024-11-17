use std::path::PathBuf;
use clap::Args;
use anyhow::{
    Context,
    anyhow,
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
        },
    }
};


type c64 = Complex<f64>;
type c32 = Complex<f32>;


#[derive(Debug, Args)]
/// Calculate model non-adiabatic coupling (NAC) for NAMD-LMI (a subset of Hefei-NAMD).
///
/// This NAC contains eigenvalues and transition dipole moment (TDM) calculated from selected WAVECAR.
/// The phonon contribution is wipped out, which means the eigenvalues and TDM stay staic over the
/// time, and the NAC (<i|d/dt|j>) vanishes.
///
/// Detailed fields of the produced file:{n}
/// - 
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

    #[arg(long, num_args(2))]
    /// Selected band range, starts from 1.
    ///
    /// Example: --brange 24 42
    brange: Vec<usize>,

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
                "x" => WavecarType::GamaHalf(Axis::X),
                "z" => WavecarType::GamaHalf(Axis::Z),
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
        let nspin   = wav.nspin;
        let ikpoint = self.ikpoint;
        let nbands  = wav.nbands;
        let nbrange = brange[1] - brange[0] + 1;

        let olaps = na::Array4::<f64>::zeros((nsw, nspin as usize, nbrange, nbrange));
        let mut eigs  = na::Array3::<f64>::zeros((nsw, nspin as usize, nbrange));
        let mut pijs  = na::Array5::<c64>::zeros((nsw, nspin as usize, 3, nbrange, nbrange));

        //wav.band_eigs.slice(na::s![])
        for isw in  0..nsw {
            eigs.slice_mut(na::s![isw, .., ..])
                .assign(&wav.band_eigs.slice(na::s![.., ikpoint-1, brange[0]-1 .. brange[1]]));
        }

        for ispin in 0 .. nspin {
            //let phi_i = wav.rekk

            //for (ii, iband) in (brange[0] - 1 .. brange[1]).enumerate() {
                //for (jj, jband) in (brange[0] - 1 .. brange[1]).enumerate().skip(ii+1) {
                    ////let pij = wav.transition_dipole_moment
                    ////pijs.slice_mut(na::s![isw, ispin])
                //}
            //}
        }

        Ok(())
    }
}
