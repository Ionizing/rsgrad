use std::path::PathBuf;

use structopt::StructOpt;

use crate::{
    vasp_parsers::wavecar::Wavecar,
    types::{
    Result,
    OptProcess,
    },
};

#[derive(Debug, StructOpt)]
pub struct Wav {
    #[structopt(long)]
    dummy: bool,
}

impl OptProcess for Wav {
    fn process(&self) -> Result<()> {
        Ok(())
    }
}
