use std::path::PathBuf;
use std::str::FromStr;

use clap::Args;
use log::info;
use anyhow::Result;

use crate::{
    OptProcess,
    pawpot::{PawPoscar, PawWavecar, PawPotcar, compute_overlap, write_overlap, write_overlap_npz},
};


#[derive(Debug, Args)]
/// Compute AE wavefunction overlaps between two WAVECARs.
///
/// Computes S_{ij}(k) = <Φ_i^a(k)|Φ_j^b(k)> including the PAW one-centre correction.
pub struct Overlap {
    #[arg(long, default_value = "WAVECAR")]
    /// First WAVECAR file path.
    wavecar1: PathBuf,

    #[arg(long)]
    /// Second WAVECAR file path (required).
    wavecar2: PathBuf,

    #[arg(long, default_value = "POSCAR")]
    /// First POSCAR file path.
    poscar1: PathBuf,

    #[arg(long)]
    /// Second POSCAR file path (required).
    poscar2: PathBuf,

    #[arg(long, default_value = "POTCAR")]
    /// POTCAR file path (shared between both structures).
    potcar: PathBuf,

    #[arg(long, num_args(1..), default_values_t = vec![1usize])]
    /// Spin indices (1-indexed). Multiple values allowed.
    ispins: Vec<usize>,

    #[arg(long, num_args(1..), default_values_t = vec![1usize])]
    /// K-point indices (1-indexed). Multiple values allowed.
    kpoints: Vec<usize>,

    #[arg(long, num_args(1..), required = true)]
    /// Initial band index specs (e.g. "1:10" "15" "1:20:2"). 1-indexed.
    ibands: Vec<String>,

    #[arg(long, num_args(1..), required = true)]
    /// Final band index specs (e.g. "1:10" "15" "1:20:2"). 1-indexed.
    jbands: Vec<String>,

    #[arg(long)]
    /// If set, skip the PAW one-centre correction (pseudo-only overlap).
    pseudo: bool,

    #[arg(short, long, default_value = "overlap.dat")]
    /// Output file path. Extension determines format: .npz → numpy, else text table.
    output: PathBuf,
}


impl OptProcess for Overlap {
    fn process(&self) -> Result<()> {
        info!("Reading POSCAR1 from {:?}", &self.poscar1);
        let poscar1 = PawPoscar::from_file(&self.poscar1)?;

        info!("Reading POSCAR2 from {:?}", &self.poscar2);
        let poscar2 = PawPoscar::from_file(&self.poscar2)?;

        info!("Reading POTCAR from {:?}", &self.potcar);
        let pawpot = PawPotcar::from_file(&self.potcar)?;

        info!("Reading WAVECAR1 from {:?}", &self.wavecar1);
        let mut wavecar1 = PawWavecar::from_file(&self.wavecar1)?;

        info!("Reading WAVECAR2 from {:?}", &self.wavecar2);
        let mut wavecar2 = PawWavecar::from_file(&self.wavecar2)?;

        let ibands = parse_band_indices(&self.ibands)?;
        let jbands = parse_band_indices(&self.jbands)?;

        info!("Computing overlaps...");
        let entries = compute_overlap(
            &poscar1,
            &poscar2,
            &pawpot,
            &mut wavecar1,
            &mut wavecar2,
            &self.ispins,
            &self.kpoints,
            &ibands,
            &jbands,
            self.pseudo,
        )?;

        let ext = self.output.extension()
            .and_then(|e| e.to_str())
            .unwrap_or("")
            .to_lowercase();

        match ext.as_str() {
            "npz" => {
                info!("Writing overlaps to {:?} as .npz ...", &self.output);
                write_overlap_npz(&self.output, &entries)?;
            }
            _ => {
                info!("Writing overlaps to {:?} as text table ...", &self.output);
                write_overlap(&self.output, &entries)?;
            }
        }

        info!("Done.");
        Ok(())
    }
}


fn parse_band_indices(specs: &[String]) -> Result<Vec<usize>> {
    let mut bands = Vec::new();
    for spec in specs {
        let parts: Vec<&str> = spec.split(':').collect();
        match parts.len() {
            1 => bands.push(usize::from_str(parts[0])?),
            2 => {
                let start = usize::from_str(parts[0])?;
                let end = usize::from_str(parts[1])?;
                bands.extend(start..=end);
            }
            3 => {
                let start = usize::from_str(parts[0])?;
                let end = usize::from_str(parts[1])?;
                let step = usize::from_str(parts[2])?;
                bands.extend((start..=end).step_by(step));
            }
            _ => anyhow::bail!("invalid band spec: {spec}"),
        }
    }
    Ok(bands)
}
