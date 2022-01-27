use std::{
    path::{Path, PathBuf},
    io::Write,
    fs,
    fmt,
};
use colored::Colorize;
use log::{
    info,
    warn,
    debug,
};
use structopt::{
    StructOpt,
    clap::AppSettings,
};
use crate::{
    traits::{
        Result,
        OptProcess
    },
    vasp_parsers::outcar::Outcar
};


#[derive(Debug, StructOpt)]
#[structopt(setting = AppSettings::ColoredHelp,
            setting = AppSettings::ColorAuto)]
/// Tracking relaxation or MD progress.
///
/// Contains the evolution of energy, maximum of Hellmann-Feynman forces, 
/// magnetic moments and time usage of each ionic step.
///
/// Hint: This command may require POSCAR for atom constraints information.
pub struct Rlx {
    #[structopt(short = "o", long = "outcar", default_value = "./OUTCAR")]
    outcar: PathBuf,

    #[structopt(short = "p", long = "poscar", default_value = "./POSCAR")]
    poscar: PathBuf,

    #[structopt(short = "e", long = "toten")]
    /// Prints TOTEN in eV
    print_energy: bool,

    #[structopt(short = "a", long = "favg")]
    /// Prints averaged total force in eV/A
    print_favg: bool,

    #[structopt(short = "x", long = "fmaxis")]
    /// Prints the axis where the strongest total force component lies on. [XYZ]
    print_fmax_axis: bool,

    #[structopt(short = "i" ,long = "fmidx")]
    /// Prints the index of ion with maximum total force load. Starts from 1
    print_fmax_index: bool,

    #[structopt(short = "v", long = "volume")]
    /// Prints lattice volume in A^3
    print_volume: bool,

    #[structopt(long = "no-fmax")]
    /// Don't print maximum total force in A^3
    no_print_fmax: bool,

    #[structopt(long = "no-totenz")]
    /// Don't print TOTEN without entropy in eV
    no_print_energyz: bool,

    #[structopt(long = "no-lgde")]
    /// Don't print Log10(delta(TOTEN without entropy))
    no_print_lgde: bool,

    #[structopt(long = "no-magmom")]
    /// Don't print total magnetic moment in muB
    no_print_magmom: bool,

    #[structopt(long = "no-nscf")]
    /// Don't print number of SCF iteration for each ionic step
    no_print_nscf: bool,

    #[structopt(long = "no-time")]
    /// Don't print time elapsed for each ionic step in minutes
    no_print_time: bool,
}


impl OptProcess for Rlx {
    fn process(&self) -> Result<()> {
        info!("Parsing file \"{:?}\" and \"{:?}\"", &self.outcar, &self.poscar);
        debug!("    OUTCAR file path = {:?}\n    POSCAR file path = {:?}",
               fs::canonicalize(&self.outcar), fs::canonicalize(&self.poscar));

        let mut outcar = Outcar::from_file(&self.outcar)?;
        if let Ok(poscar) = vasp_poscar::Poscar::from_path(&self.poscar) {
            if let Some(constraints) = poscar.into_raw().dynamics {
                outcar.set_constraints(constraints);
            }
        } else {
            warn!("Reading connstraints from POSCAR file \"{:?}\" failed", &self.poscar);
        }

        Ok(())
    }
}