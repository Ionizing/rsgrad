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
        let fname = PathBuf::from("WAVECAR");
        let wav = Wavecar::from_file(&fname)?;

        println!("{}", &wav);
        println!("{:#}", &wav);

        Ok(())
    }
}
