use std::path::PathBuf;

use clap::Args;
use log::info;

use crate::{
    types::Result,
    OptProcess,
    pawpot::{PawPoscar, PawWavecar, PawPotcar, write_normalcar, write_cproj_npz},
};


#[derive(Debug, Args)]
/// Write NormalCAR (PAW projector coefficients) from POSCAR + POTCAR + WAVECAR.
///
/// The NormalCAR stores the PAW projector coefficients β = ⟨p̃_i|ψ̃⟩ for every
/// spin/k-point/band, useful for spinorb and other post-processing tools.
pub struct Normalcar {
    #[arg(long, default_value = "WAVECAR")]
    /// WAVECAR file path.
    wavecar: PathBuf,

    #[arg(long, default_value = "POSCAR")]
    /// POSCAR file path.
    poscar: PathBuf,

    #[arg(long, default_value = "POTCAR")]
    /// POTCAR file path.
    potcar: PathBuf,

    #[arg(short, long, default_value = "NormalCAR")]
    /// Output file path. If extension is .npz, write numpy npz format; otherwise write binary NormalCAR.
    output: PathBuf,
}


impl OptProcess for Normalcar {
    fn process(&self) -> Result<()> {
        info!("Reading POSCAR from {:?}", &self.poscar);
        let poscar = PawPoscar::from_file(&self.poscar)?;

        info!("Reading POTCAR from {:?}", &self.potcar);
        let pawpot = PawPotcar::from_file(&self.potcar)?;

        info!("Reading WAVECAR from {:?}", &self.wavecar);
        let mut wavecar = PawWavecar::from_file(&self.wavecar)?;

        let ext = self.output.extension()
            .and_then(|e| e.to_str())
            .unwrap_or("")
            .to_lowercase();

        match ext.as_str() {
            "npz" => {
                info!("Writing cproj to {:?} as .npz ...", &self.output);
                write_cproj_npz(&self.output, &poscar, &pawpot, &mut wavecar)?;
            }
            _ => {
                info!("Writing NormalCAR to {:?} ...", &self.output);
                write_normalcar(&self.output, &poscar, &pawpot, &mut wavecar)?;
            }
        }

        info!("Done.");
        Ok(())
    }
}
