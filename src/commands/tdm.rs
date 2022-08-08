use std::path::PathBuf;
use clap::{
    Parser,
    AppSettings,
};

use crate::{
    types::Result,
    OptProcess,
    vasp_parsers::wavecar::Wavecar
};

#[derive(Debug, Parser)]
#[clap(setting = AppSettings::ColoredHelp,
       setting = AppSettings::ColorAuto)]
/// Calculate Transition Dipole Moment (TDM) between given bands.
///
/// Note: This command can only calculate the TDM between bands in
/// same k-point. The inter-kpoint transition is not supported yet.
/// Also, this commands calculates the TDM in reciprocal space by{n}
///
/// tdm_{i->j} = <phi_j|e*r|phi_i> = i*ħ/(ΔE*m)*<phi_j|p|phi_i>
pub struct Tdm {
    #[clap(short, long, default_value = "./WAVECAR")]
    /// WAVECAR file path.
    wavecar: PathBuf,

    #[clap(short = 's', long, default_value = "1", possible_values = &["1", "2"])]
    /// Spin index, 1 for up, 2 for down.
    ispin: usize,

    #[clap(short = 'k', long, default_value = "1")]
    /// K-point index, starts from 1.
    ikpoint: usize,

    #[clap(short = 'i', long, min_values = 1, required = true)]
    /// Initial band indices, start from 1.
    ibands: Vec<usize>,

    #[clap(short = 'j', long, min_values = 1, required = true)]
    /// Final band indices, starts from 1.
    jbands: Vec<usize>,

    #[clap(long, default_value = "0.05")]
    /// Smearing width, in eV.
    sigma: f64,

    #[clap(long, default_value = "tdm_peaks.txt")]
    /// Write the TDM peaks to raw txt file.
    peakout: PathBuf,

    #[clap(long, default_value = "tdm_smeared.txt")]
    /// Write the summed and smeared TDM to raw txt file.
    txtout: PathBuf,

    #[clap(long, default_value = "tdm_smeared.html")]
    /// Write the plot of TDM to html file.
    htmlout: PathBuf,
}


impl OptProcess for Tdm {
    fn process(&self) -> Result<()> {
        Ok(())
    }
}
