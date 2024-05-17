use std::{
    path::PathBuf,
    fs,
};
use rayon::prelude::*;
use log::{
    info,
    warn,
    debug,
};
use clap::Args;
use crate::{
    Result,
    index_transform,
    OptProcess,
    Outcar,
    Poscar,
    Vibrations,
    vasp_parsers::outcar::PrintAllVibFreqs,
};


#[derive(Debug, Args)]
#[command(allow_negative_numbers = true)]
/// Tracking vibration information.
///
/// For systems enabled vibration mode calculation, this command can extract
/// phonon eigenvalues and phonon eigenvectors at Gamma point.
pub struct Vib {
    #[arg(short = 'l', long)]
    /// Shows vibration modes in brief
    list: bool,

    #[arg(default_value = "./OUTCAR")]
    /// Specify the input OUTCAR file
    outcar: PathBuf,

    #[arg(short = 'p', long, default_value = "./POSCAR")]
    /// Specify the input POSCAR file, the consntraints info is needed
    poscar: PathBuf,

    #[arg(short = 'x', long)]
    /// Saves each selected modes to XSF file
    save_as_xsfs: bool,

    #[arg(short = 'i', long, num_args(0..))]
    /// Selects the mode indices to operate.
    ///
    /// Step indices start from '1', if '0' is given, all the structures will be selected.
    /// Step indices can be negative, where negative index means counting reversely.
    /// E.g. "-i -2 -1 1 2 3" means selecting the last two and first three steps.
    select_indices: Option<Vec<i32>>,

    #[arg(long, default_value = ".")]
    /// Define where the files would be saved
    save_in: PathBuf,

    #[arg(short = 'm', long)]
    /// Modulate the ground-state POSCAR with respect to a certern vibration frequencies.
    modulate: bool,

    #[arg(short = 'a', long, default_value_t = 0.1)]
    /// Modulation amplitude coefficient, to avoid precision issue, abs(amplitude) >= 0.01 should
    /// be satisfied.
    amplitude: f64,

    #[arg(short = 'c', long)]
    /// Catesian coordinate is used when writting POSCAR. Fractional coordinate is used by default.
    cartesian: bool,
}


impl OptProcess for Vib {
    fn process(&self) -> Result<()> {
        info!("Parsing file {:?}", &self.outcar);
        debug!("    OUTCAR file path = {:?}", fs::canonicalize(&self.outcar));

        let outcar = Outcar::from_file(&self.outcar)?;
        let vibs = Vibrations::from(outcar);

        let constraints = if let Ok(poscar) = Poscar::from_file(&self.poscar) {
            poscar.constraints
        } else {
            warn!("Reading constraints from POSCAR file {:?} failed", &self.poscar);
            None
        };

        let len = vibs.modes.len();
        let inds: Vec<usize> = index_transform(self.select_indices.clone().unwrap_or_default(), len);

        if self.list {
            let paf: PrintAllVibFreqs = vibs.clone().into();
            print!("{}", paf);
        }

        if self.save_as_xsfs {
            inds.par_iter()
                .map(|i| {
                    vibs.save_as_xsf(*i, &self.save_in)?;
                    Ok(())
                })
                .collect::<Result<()>>()?;
        }

        if self.modulate {
            inds.par_iter()
                .map(|i| {
                    vibs.modulate(*i, self.amplitude, !self.cartesian, constraints.clone(), &self.save_in)?;
                    Ok(())
                })
                .collect::<Result<()>>()?;
        }

        Ok(())
    }
}
