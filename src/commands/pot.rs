use std::path::PathBuf;

use clap::Args;
use log::{
    info,
    warn,
};
use anyhow::Context;

use crate::{
    Poscar,
    Potcar,
    Settings,
    Result,
    OptProcess,
    vasp_parsers::potcar::FunctionalType,
};


#[derive(Debug, Args)]
/// Generate the POTCAR according to POSCAR
pub struct Pot {
    #[arg(long, short)]
    /// Specify the configuration file, if left blank, rsgrad will read `.rsgrad.toml` at
    /// your home dir.
    config: Option<PathBuf>,

    #[arg(long, short, default_value = "./POSCAR")]
    /// Specify the POSCAR file
    ///
    /// Note: the elements symbols' line should involve the specified valence configuration
    /// if you need. For example, if you need a `K_sv` POTCAR for `K`, just write `K_sv` in
    /// the POSCAR
    poscar: PathBuf,

    #[arg(default_value = "PAW_PBE")]
    /// Specify the functional type, now only "PAW_PBE"(or "paw_pbe") and "PAW_LDA"(or "paw_lda") are available.
    functional: FunctionalType,

    #[arg(long, short, default_value = "./")]
    /// Specify where the `POTCAR` would be written
    save_in: PathBuf,
}


impl OptProcess for Pot {
    fn process(&self) -> Result<()> {
        let settings = if let Some(path) = self.config.as_ref() {
            Settings::from_file(path)?
        } else {
            Settings::from_default()?
        };

        info!("Reading POSCAR file {:?} ...", &self.poscar);
        let pos = Poscar::from_file(&self.poscar)?;

        let (symbols, specified_types) = {
            let titels = &pos.ion_types;
            let mut symbols = Vec::<String>::new();
            let mut specified_types = Vec::<String>::new();
            
            for titel in titels {
                let mut it = titel.splitn(2, '_');
                let symbol = it.next()
                    .context(format!("No chemical symbol found in {}", titel))?
                    .to_string();
                let specified_type = it.next()
                    .map(|x| "_".to_string() + x )
                    .unwrap_or_default();

                symbols.push(symbol);
                specified_types.push(specified_type);
            }

            (symbols, specified_types)
        };

        let pot = Potcar::from_config(&symbols,
                                      &self.functional,
                                      &specified_types,
                                      &settings.functional_path)?;

        let fname = self.save_in.with_file_name("POTCAR");
        if fname.is_file() {
            warn!("Found `POTCAR` in specified dir, renaming it to `POTCAR.bak`");
            let renamed_fname = fname.with_extension("bak");
            std::fs::rename(&fname, renamed_fname)?;
        }

        info!("Writing POTCAR to {:?} ...", &fname);
        pot.to_file(&fname)?;

        Ok(())
    }
}
