use std::{
    io::prelude::*,
    path::Path,
    fs::read_to_string,
};

use anyhow::{
    Context,
    Result,
    bail,
};
use flate2::read::GzDecoder;
use log::info;

use crate::settings::FunctionalPath;


#[allow(non_camel_case_types)]
#[derive(Clone, Copy, Debug)]
pub enum FunctionalType {
    PAW_PBE,
    PAW_LDA,
}


impl std::fmt::Display for FunctionalType {
    fn fmt(&self, f: &mut std::fmt::Formatter<'_>) -> std::fmt::Result {
        let s = match self {
            FunctionalType::PAW_PBE => "PAW_PBE",
            FunctionalType::PAW_LDA => "PAW_LDA",
        };
        write!(f, "{}", s)
    }
}


impl std::str::FromStr for FunctionalType {
    type Err = anyhow::Error;

    fn from_str(s: &str) -> Result<Self> {
        let ret = match s {
            "PAW_PBE" | "paw_pbe" => Self::PAW_PBE,
            "PAW_LDA" | "paw_lda" => Self::PAW_LDA,
            _ => {
                bail!(r#"``{}` cannot be converted into FunctionalType.
Available functionals are `PAW_PBE`(or `paw_pbe`) and `PAW_LDA`(or `paw_lda`)."#, s);
            }
        };
        Ok(ret)
    }
}


#[derive(Clone, Debug)]
pub struct AtomicPotcar {
    pub symbol: String,                 // Element symbol, H, He, Li, Be, B, C ...
    pub functional: FunctionalType,     // Functional type, LDA, PBE
    pub specific_type: String,          // Valence annotation '_sv', '_GW', '_AE' ...
    pub content: String,                // Raw content of single element POTCAR
}


impl AtomicPotcar {
    pub fn from_config(symbol: &str, 
                       functional: &FunctionalType, 
                       specific_type: &str,
                       prefix: &FunctionalPath) -> Result<Self> {
        let titel = symbol.to_string() + specific_type;

        //let titel = prefix.aliases.as_ref()
            //.and_then(|tok| { tok.get(&titel).cloned() })
            //.unwrap_or(titel);
        let titel = if let Some(alias) = prefix.aliases.as_ref().and_then(|tok| tok.get(&titel)) {
            info!("Found functional alias: \"{}\" -> \"{}\".", titel, alias);
            alias.clone()
        } else {
            titel
        };

        let path = {
            let mut _ret = match functional {
                FunctionalType::PAW_PBE => prefix.paw_pbe.to_path_buf(),
                FunctionalType::PAW_LDA => prefix.paw_lda.to_path_buf(),
            };
            _ret.push(&titel);
            _ret.push("POTCAR");

            _ret.canonicalize().unwrap()
        };

        info!("Reading POTCAR from {:?}", &path);

        let content = if path.is_file() {
            read_to_string(&path)?
        } else {
            let fname = vec![
                path.with_extension(".z"),
                path.with_extension(".Z"),
                path.with_extension(".gz"),
            ].into_iter().find(|p| p.is_file())
                .context(format!("No suitable POTCAR found for element {}", symbol))?;

            let bytes = std::fs::read(fname)?;
            let mut gz = GzDecoder::new(&bytes[..]);
            let mut s = String::new();
            gz.read_to_string(&mut s)?;

            s
        };

        Ok(
            Self {
                symbol: symbol.to_string(),
                functional: *functional,
                specific_type: specific_type.to_string(),
                content
            }
            )
    }
}


pub struct Potcar {
    pub inner: Vec<AtomicPotcar>,
}


impl Potcar {
    pub fn from_config(symbols: &[String],
                       functional: &FunctionalType,
                       specific_types: &[String],
                       prefix: &FunctionalPath) -> Result<Self> {
        let mut inner = Vec::<AtomicPotcar>::new();

        for (sym, spec) in symbols.iter().zip(specific_types.iter()) {
            inner.push(AtomicPotcar::from_config(sym, functional, spec, prefix)?);
        }

        Ok(Self { inner })
    }


    /// Concatenate all the AtomicPotcar.content into a String
    pub fn to_txt(&self) -> String {
        self.inner.iter()
            .fold(String::new(), |mut accum, item| {
                accum.push_str(&item.content);
                accum
            })
    }


    pub fn to_file(&self, path: &(impl AsRef<Path> + ?Sized)) -> Result<()> {
        std::fs::write(path, self.to_txt())?;
        Ok(())
    }
}


#[cfg(test)]
mod tests {
    use super::*;
    use crate::settings::Settings;

    #[test]
    #[ignore]
    fn test_atomic_potcar() {
        let symbol = "K";
        let functional = FunctionalType::PAW_PBE;
        let specific_type = "_sv";
        let prefix = &Settings::from_default().unwrap().functional_path;

        let potcar = AtomicPotcar::from_config(symbol, &functional, specific_type, prefix)
            .unwrap()
            .content;
        print!("{}", potcar);
    }

    #[test]
    #[ignore]
    fn test_potcar_from_config() {
        let symbols = vec!["Ag".to_owned(), "K".to_owned(), "N".to_owned()];
        let functional = FunctionalType::PAW_PBE;
        let specific_types = vec!["".to_owned(), "_sv".to_owned(), "".to_owned()];
        let prefix = &Settings::from_default().unwrap().functional_path;

        let potcar = Potcar::from_config(&symbols, &functional, &specific_types, prefix).unwrap().to_txt();
        print!("{}", potcar);
    }
}
