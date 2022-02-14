use std::{
    io::prelude::*,
    fs::read_to_string,
};

use anyhow::{
    Context,
    Result,
};
use flate2::read::GzDecoder;

use crate::settings::FunctionalPath;


#[allow(non_camel_case_types)]
#[derive(Clone, Copy, Debug)]
pub enum FunctionalType {
    PAW_PBE,
    PAW_LDA,
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
        let path = {
            let mut _ret = match functional {
                FunctionalType::PAW_PBE => prefix.paw_pbe.to_path_buf(),
                FunctionalType::PAW_LDA => prefix.paw_lda.to_path_buf(),
            };
            _ret.push(&titel);
            _ret.push("POTCAR");

            _ret.canonicalize().unwrap()
        };

        println!("{:?}", &path);

        let content = if path.is_file() {
            read_to_string(&path)?
        } else {
            let fname = vec![
                path.with_extension(".z"),
                path.with_extension(".Z"),
                path.with_extension(".gz"),
            ].into_iter().find(|p| p.is_file())
                .context(format!("No suitable POTCAR found for element {}", symbol))?;

            let bytes = std::fs::read(&fname)?;
            let mut gz = GzDecoder::new(&bytes[..]);
            let mut s = String::new();
            gz.read_to_string(&mut s)?;

            s
        };

        Ok(
            Self {
                symbol: symbol.to_string(),
                functional: functional.clone(),
                specific_type: specific_type.to_string(),
                content
            }
            )
    }
}


pub struct Potcar {
    pub inner: Vec<AtomicPotcar>,
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
        print!("{}", AtomicPotcar::from_config(symbol, &functional, specific_type, &prefix).unwrap()
               .content);
    }
}
