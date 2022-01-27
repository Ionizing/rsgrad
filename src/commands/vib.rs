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
    vasp_parsers::outcar
};


#[derive(Debug, StructOpt)]
#[structopt(setting = AppSettings::ColoredHelp,
            setting = AppSettings::ColorAuto)]
/// Tracking vibration information.
///
/// For systems enabled vibration mode calculation, this command can extract
/// phonon eigenvalues and phonon eigenvectors at Gamma point.
pub struct Vib {
    #[structopt(short, long)]
    /// Shows vibration modes in brief
    list: bool,

    #[structopt(short = "x", long)]
    /// Saves each selected modes to XSF file
    save_as_xsfs: bool,

    #[structopt(short = "i", long)]
    /// Selects the indices to operate.
    ///
    /// Step indices start from '1', if '0' is given, all the structures will be selected.
    /// Step indices can be negative, where negative index means counting reversely.
    /// E.g. "--save-as-poscars -2 -1 1 2 3" means saving the last two and first three
    /// steps.
    select_indices: Option<Vec<i32>>,

    #[structopt(long, default_value = ".")]
    /// Define where the files would be saved
    save_in: PathBuf,
}


impl OptProcess for Vib {
    fn process(&self) -> Result<()> {
        Ok(())
    }
}
