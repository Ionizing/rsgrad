use std::path::PathBuf;

use clap::Args;
use log::info;

use crate::{
    types::Result,
    OptProcess,
    pawpot::{PawPoscar, PawWavecar, PawPotcar, VaspAeWfc},
};


#[derive(Debug, Args)]
/// Reconstruct the all-electron (AE) wavefunction density from a WAVECAR.
///
/// Uses the PAW transformation to reconstruct the full AE wavefunction from
/// the pseudo-wavefunction stored in the WAVECAR.
pub struct Aewfc {
    #[arg(long, default_value = "WAVECAR")]
    /// WAVECAR file path.
    wavecar: PathBuf,

    #[arg(long, default_value = "POSCAR")]
    /// POSCAR file path.
    poscar: PathBuf,

    #[arg(long, default_value = "POTCAR")]
    /// POTCAR file path.
    potcar: PathBuf,

    #[arg(long, default_value_t = 1)]
    /// Spin index, 1-indexed.
    ispin: usize,

    #[arg(long, default_value_t = 1)]
    /// K-point index, 1-indexed.
    ikpt: usize,

    #[arg(long, default_value_t = 1)]
    /// Band index, 1-indexed.
    iband: usize,

    #[arg(long, default_value_t = -2.0)]
    /// AE energy cutoff in eV. Negative values mean |aecut| * pscut.
    aecut: f64,

    #[arg(short, long, default_value = "AE_density.vasp")]
    /// Output file path. Extension determines format: .npy (wfc array), .npz (wfc+density), else VASP CHGCAR.
    output: PathBuf,
}


impl OptProcess for Aewfc {
    fn process(&self) -> Result<()> {
        info!("Reading POSCAR from {:?}", &self.poscar);
        let poscar = PawPoscar::from_file(&self.poscar)?;

        info!("Reading POTCAR from {:?}", &self.potcar);
        let pawpot = PawPotcar::from_file(&self.potcar)?;

        info!("Reading WAVECAR from {:?}", &self.wavecar);
        let wavecar = PawWavecar::from_file(&self.wavecar)?;

        info!("Building AE-WFC reconstruction object for k-point {}...", self.ikpt);
        let mut aewfc = VaspAeWfc::new(wavecar, &poscar, &pawpot, self.ikpt, self.aecut)?;

        let ext = self.output.extension()
            .and_then(|e| e.to_str())
            .unwrap_or("")
            .to_lowercase();

        match ext.as_str() {
            "npy" => {
                info!("Writing AE wavefunction to {:?} as .npy ...", &self.output);
                aewfc.write_ae_wfc_npy(self.ispin, self.iband, &self.output)?;
            }
            "npz" => {
                info!("Writing AE wavefunction + density to {:?} as .npz ...", &self.output);
                aewfc.write_ae_wfc_npz(self.ispin, self.iband, &self.output)?;
            }
            _ => {
                info!("Writing AE density to {:?} as VASP CHGCAR format ...", &self.output);
                aewfc.write_ae_density(&poscar, self.ispin, self.iband, &self.output)?;
            }
        }

        info!("Done.");
        Ok(())
    }
}
