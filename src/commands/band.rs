use std::{
    fs,
    path::PathBuf,
};

use structopt::{
    StructOpt,
    clap::AppSettings,
};
use log::{
    info,
    debug,
};

use crate::{
    Result,
    OptProcess,
    Procar,
    vasp_parsers::outcar::GetEFermi,
    types::Vector,
    commands::common::write_array_to_txt,
};

struct Configuration {

}


pub struct Band {
    config: Configuration,

    gen_template: bool,

    procar: PathBuf,

    outcar: PathBuf,

    txtout: PathBuf,

    htmlout: PathBuf,
}
