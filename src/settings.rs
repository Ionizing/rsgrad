use std::path::{
    Path,
    PathBuf,
};
use anyhow::{
    Result,
    Context,
    bail,
};
use log::{
    info,
    debug,
};
use config::Config;
use directories::BaseDirs;
use serde::{
    Deserialize,
    Serialize,
};
use colored::Colorize;


#[derive(Debug, Deserialize, Serialize, PartialEq, Eq)]
pub struct Settings {
    #[serde(rename(serialize   = "functional-path",
                   deserialize = "functional-path"))]
    functional_path: FunctionalPath,
}


#[derive(Debug, Deserialize, Serialize, PartialEq, Eq)]
pub struct FunctionalPath {
    #[serde(rename(serialize   = "PAW_PBE",
                   deserialize = "PAW_PBE"))]
    pub paw_pbe: PathBuf,

    #[serde(rename(serialize   = "PAW_LDA",
                   deserialize = "PAW_LDA"))]
    pub paw_lda: PathBuf,
}


impl Settings {
    pub fn from_file(path: &(impl AsRef<Path> + ?Sized)) -> Result<Self> {
        info!("Reading rsgrad settings from {:?} ...", path.as_ref());
        Self::check_file_availability(path)?;
        
        let settings = Config::builder()
            .add_source(config::File::new(path.as_ref().to_str().unwrap(), config::FileFormat::Toml))
            .build()?
            .try_deserialize::<Settings>()?;
        debug!("Configuration file {:?} content: {:#?}", path.as_ref(), settings);
        settings.check_availability()?;

        Ok(settings)
    }

    pub fn from_config() -> Result<Self> {
        let mut path: PathBuf = BaseDirs::new()
            .context("Home directory not found.")?
            .home_dir().to_path_buf();
        path.push(".rsgrad.toml");

        if !path.is_file() {
            let help_conf = r#"[functional-path]
PAW_PBE = "/public/apps/vasp/potpaw_PBE.54"
PAW_LDA = "/public/apps/vasp/potpaw_LDA.54""#.bright_yellow();

            let example_conf = r#"[functional-path]
PAW_PBE = "<path of PAW_PBE>"
PAW_LDA = "<path of PAW_LDA>""#.bright_yellow();

            bail!(r#"rsgrad configuration file {:?} is not a regular file or doesn't exist.
Consider create that file with similar content in the following:

{}

Please replace {} with actual path of corresponding PP's directory, for example:

{}"#, path, example_conf, "<path of ...>".bright_yellow(), help_conf);
        }

        Self::from_file(&path)
    }

    fn check_availability(&self) -> Result<()> {
        info!("Checking rsgrad setting availability ...");

        Self::check_dir_availability(&self.functional_path.paw_pbe)?;
        Self::check_dir_availability(&self.functional_path.paw_lda)?;

        Ok(())
    }

    fn check_dir_availability(dir: &(impl AsRef<Path> + ?Sized)) -> Result<()> {
        if !dir.as_ref().is_dir() {
            bail!("Directory {:?} not available. It should be a regular directory.", dir.as_ref())
        } else {
            Ok(())  
        }
    }

    fn check_file_availability(path: &(impl AsRef<Path> + ?Sized)) -> Result<()> {
        if !path.as_ref().is_file() {
            bail!("File {:?} not available. It should be a regular file.", path.as_ref())
        } else {
            Ok(())
        }
    }
}


#[cfg(test)]
mod tests {
    use super::*;
    use toml;
    
    #[test]
    fn test_serialize() {
        let settings = Settings {
            functional_path: FunctionalPath {
                paw_pbe: PathBuf::from("~/apps/vasp/potpaw_PBE.54"),
                paw_lda: PathBuf::from("~/apps/vasp/potpaw_LDA.54"),
            }
        };

        let txt = r#"[functional-path]
PAW_PBE = "~/apps/vasp/potpaw_PBE.54"
PAW_LDA = "~/apps/vasp/potpaw_LDA.54"
"#;
        assert_eq!(toml::to_string(&settings).unwrap(), txt);
    }

    #[test]
    fn test_deserialize() {
        let settings_expected = Settings {
            functional_path: FunctionalPath {
                paw_pbe: PathBuf::from("~/apps/vasp/potpaw_PBE.54"),
                paw_lda: PathBuf::from("~/apps/vasp/potpaw_LDA.54"),
            }
        };

        let txt = r#"[functional-path]
PAW_PBE = "~/apps/vasp/potpaw_PBE.54"
PAW_LDA = "~/apps/vasp/potpaw_LDA.54"
"#;

        let parsed: Settings = toml::from_str(txt).unwrap();
        assert_eq!(parsed, settings_expected);
    }

    #[test]
    #[ignore]
    fn test_from_config() -> Result<()> {
        Settings::from_config()?;
        Ok(())
    }
}
