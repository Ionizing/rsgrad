use std::path::PathBuf;

use structopt::{
    StructOpt,
    clap::AppSettings,
};
use anyhow::{
    Error,
    bail,
};

use crate::{
    Result,
    OptProcess,
    Procar
};


#[derive(Debug)]
struct Selections(Vec<(Vec<i32>, String)>);


impl std::str::FromStr for Selections {
    type Err = anyhow::Error;
    fn from_str(s: &str) -> Result<Self> {

        bail!("");
    }
}



#[derive(Debug, StructOpt)]
#[structopt(setting = AppSettings::ColoredHelp,
            setting = AppSettings::ColorAuto)]
pub struct Dos {
    config: Option<PathBuf>,

    outcar: PathBuf,
    procar: PathBuf,
    ikpoints: Vec<i32>,

    projections: Selections,

    txtout: PathBuf,
    htmlout: PathBuf,
}


impl OptProcess for Dos {
    fn process(&self) -> Result<()> {
        Ok(())
    }
}
