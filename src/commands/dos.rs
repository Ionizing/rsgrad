use std::path::PathBuf;

use structopt::{
    StructOpt,
    clap::AppSettings,
};
use serde;
use anyhow::{
    Error,
    bail,
};

use crate::{
    Result,
    OptProcess,
    Procar
};


struct Selection {
    ispins:   Vec<i32>,
    ikpoints: Vec<i32>,
    iatoms:   Vec<i32>,
    iorbits:  Vec<i32>,
    label:    String,
}


#[derive(Debug, StructOpt, Clone)]
#[structopt(setting = AppSettings::ColoredHelp,
            setting = AppSettings::ColorAuto)]
pub struct Dos {
    #[structopt(short, long)]
    config: Option<PathBuf>,

    #[structopt(long, default_value = "./OUTCAR")]
    outcar: PathBuf,

    #[structopt(long, default_value = "./PROCAR")]
    procar: PathBuf,

    #[structopt(long, default_value = "0")]
    ispins: Vec<i32>,

    #[structopt(short, long)]
    txtout: PathBuf,

    #[structopt(short, long)]
    htmlout: PathBuf,
}


impl OptProcess for Dos {
    fn process(&self) -> Result<()> {
        Ok(())
    }
}
