use std::path::PathBuf;
use clap::Args;
use log::info;
use anyhow::{
    anyhow,
    Context,
};
use rayon;
use crate::{
    types::Result,
    ChargeDensity,
    ChargeType,
    OptProcess,
};


#[derive(Debug, Args)]
/// Calculate charge density difference. 
///
/// The operation is performed by `chgdiff = chgcar_ab - (chgcar_a + chgcar_b)`.
pub struct Chgdiff {
    /// The CHGCAR of A+B system
    chgcar_ab: PathBuf,

    /// The CHGCAR of A system
    chgcar_a: PathBuf,

    /// The CHGCAR of B system
    chgcar_b: PathBuf,

    #[arg(short, long, default_value = "CHGDIFF.vasp")]
    /// The output charge density difference file path
    output: PathBuf,
}


impl OptProcess for Chgdiff {
    fn process(&self) -> Result<()> {
        let mut chgcar_ab: Result<ChargeDensity> = Err(anyhow!(""));
        let mut chgcar_a: Result<ChargeDensity>  = Err(anyhow!(""));
        let mut chgcar_b: Result<ChargeDensity>  = Err(anyhow!(""));

        rayon::scope(|s| {
            s.spawn(|_| {
                info!("Reading charge density from {:?}", self.chgcar_ab);
                chgcar_ab = ChargeDensity::from_file(&self.chgcar_ab, ChargeType::Chgcar)
                    .context(format!("Failed to read charge density from {:?}", self.chgcar_ab))
            });
            s.spawn(|_| {
                info!("Reading charge density from {:?}", self.chgcar_a);
                chgcar_a = ChargeDensity::from_file(&self.chgcar_a, ChargeType::Chgcar)
                    .context(format!("Failed to read charge density from {:?}", self.chgcar_a))
            });
            s.spawn(|_| {
                info!("Reading charge density from {:?}", self.chgcar_b);
                chgcar_b = ChargeDensity::from_file(&self.chgcar_b, ChargeType::Chgcar)
                    .context(format!("Failed to read charge density from {:?}", self.chgcar_b))
            });
        });

        let chgcar_ab = chgcar_ab?;
        let chgcar_a = chgcar_a?;
        let chgcar_b = chgcar_b?;

        info!("Calculating charge density difference by `CHGDIFF = {:?} - ({:?} + {:?})`", 
              self.chgcar_ab, self.chgcar_a, self.chgcar_b);

        let chgdiff = (chgcar_ab - (chgcar_a + chgcar_b)?)?;

        info!("Writing charge difference to {:?}", self.output);

        chgdiff.to_file(&self.output)?;

        Ok(())
    }
}
