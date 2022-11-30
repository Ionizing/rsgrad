use std::path::PathBuf;

use clap::{
    Parser,
    AppSettings,
};
use log::{
    info,
    warn,
};
use ndarray::{
    Array2,
    Array3,
};
use anyhow::Context;

use crate::{
    vasp_parsers::{
        wavecar::Wavecar,
        procar::Procar,
    },
    types::{
        Result,
        OptProcess,
    },
};


#[derive(Debug, Parser)]
#[clap(setting = AppSettings::ColoredHelp,
       setting = AppSettings::ColorAuto)]
/// Find band gap and print positions of VBM and CBM
pub struct Gap {
    #[clap(long, short = 'w', default_value = "WAVECAR")]
    /// WAVECAR file name
    wavecar: PathBuf,

    #[clap(long, short = 'p', default_value = "PROCAR")]
    /// PROCAR file name
    procar: PathBuf,
}


impl Gap {
    fn bands_from_procar(procar: &Procar) -> (Array3<f64>, Array3<f64>, Array2<f64>) {
        let eigs = procar.pdos.eigvals.clone();
        let occs = procar.pdos.occupations.clone();
        let kvec = procar.kpoints.kpointlist.clone();

        (eigs, occs, kvec)
    }

    fn bands_from_wavecar(wavecar: &Wavecar) -> (Array3<f64>, Array3<f64>, Array2<f64>) {
        let eigs = wavecar.band_eigs.clone();
        let occs = wavecar.band_fweights.clone();
        let kvec = wavecar.kvecs.clone();

        (eigs, occs, kvec)
    }

    fn get_bands_kpoints(&self) -> Result<(Array3<f64>, Array3<f64>, Array2<f64>)> {
        Procar::from_file(&self.procar).and_then(|v| Ok(Self::bands_from_procar(&v)))
            .or(Wavecar::from_file(&self.wavecar).and_then(|v| Ok(Self::bands_from_wavecar(&v))))
            .with_context(|| format!("Neither WAVECAR nor PROCAR is accessible, please specify a valid WAVECAR or PROCAR"))
    }
}


impl OptProcess for Gap {
    fn process(&self) -> Result<()> {
        let (eigs, occs, kvec) = self.get_bands_kpoints()?;

        // find cbm
        Ok(())
    }
}
