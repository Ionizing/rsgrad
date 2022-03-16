use std::{
    fs,
    path::PathBuf,
    collections::HashMap,
};

use structopt::{
    StructOpt,
    clap::AppSettings,
};
use log::{
    info,
    warn,
};
use serde::{
    Serialize,
    Deserialize,
    Deserializer,
};
use anyhow::{
    anyhow,
    Context,
    Error,
    bail,
};
use toml;

use crate::{
    Result,
    OptProcess,
    Procar,
    Outcar,
    vasp_parsers::outcar::GetEFermi,
    types::range_parse,
    types::index_transform,
};


#[derive(Clone, Serialize, Deserialize)]
struct RawSelection {
    spins:      Option<String>,
    kpoints:    Option<String>,
    atoms:      Option<String>,
    orbits:     Option<String>,
    factor:     Option<f64>,
}


#[derive(Clone, Debug)]
struct Selection {
    label:      String,
    ispins:     Vec<usize>,
    ikpoints:   Vec<usize>,
    iatoms:     Vec<usize>,
    iorbits:    Vec<usize>,
    factor:     f64,
}


fn rawsel_to_sel(r: HashMap<String, RawSelection>, 
                 nlm: &[String], 
                 nions: usize, 
                 nkpoints: usize,
                 nspin: usize,
                 is_ncl: bool) -> Result<Vec<Selection>> {

    let mut sel_vec = vec![];

    for (label, val) in r.into_iter() {

        let iorbits = if let Some(orbits) = val.orbits {
            orbits.split_whitespace()
                .map(|x| {
                    nlm.iter().position(|x2| x2 == &x)
                        .context(format!("Selected orbit {:?} not available in {:?}", x, &nlm))
                })
                .collect::<Result<Vec<_>>>()?
        } else {
            (0 .. nlm.len()).collect::<Vec<usize>>()
        };

        let ispins = if let Some(spins) = val.spins {
            if is_ncl {
                spins.split_whitespace()
                    .map(|x| match x.to_ascii_lowercase().as_ref() {
                        "x"         => Ok(1usize),
                        "y"         => Ok(2usize),
                        "z"         => Ok(3usize),
                        "t" | "tot" => Ok(0usize),
                        _ =>
bail!("Invalid spin component selected: `{}`, available components are `x`, `y`, `z` and `tot`", x)
                    }).collect::<Result<Vec<_>>>()?
            } else if nspin == 2 {
                spins.split_whitespace()
                    .map(|x| match x.to_ascii_lowercase().as_ref() {
                        "u" | "up"           => Ok(0usize),
                        "d" | "dn" | "down"  => Ok(1usize),
                        _ =>
bail!("Invalid spin component selected: `{}`, available components are `up` and `down`", x)
                    }).collect::<Result<Vec<_>>>()?
            } else {
                warn!("[DOS]: This system is not spin-polarized, only one spin component is available, selected by default.");
                vec![0usize]
            }
        } else {
            if is_ncl {
                vec![0usize]  // 'tot' part
            } else if nspin == 2 {
                vec![0usize, 1usize]  // spin up and spin down
            } else {
                vec![0usize]  // only one spin available
            }
        };

        let iatoms = if let Some(atoms) = val.atoms {
            let mut ret = atoms.split_whitespace()
                .map(|x| range_parse(x))
                .collect::<Result<Vec<Vec<i32>>>>()?
                .into_iter()
                .map(|x| index_transform(x, nions).into_iter())
                .flatten()
                .map(|x| x - 1)
                .collect::<Vec<usize>>();
            ret.sort();
            ret.dedup();

            if ret.iter().any(|x| x >= &nions) {
                bail!("[DOS]: Some selected atom index greater than NIONS");
            }

            ret
        } else {  // All atoms selected if left blank
            (0 .. nions).collect::<Vec<usize>>()
        };

        let ikpoints = if let Some(kpoints) = val.kpoints {
            let mut ret = kpoints.split_whitespace()
                .map(|x| range_parse(x))
                .collect::<Result<Vec<Vec<i32>>>>()?
                .into_iter()
                .map(|x| index_transform(x, nkpoints).into_iter())
                .flatten()
                .map(|x| x - 1)
                .collect::<Vec<usize>>();
            ret.sort();
            ret.dedup();

            if ret.iter().any(|x| x >= &nkpoints) {
                bail!("[DOS]: Some selected kpoint index greater than NKPTS");
            }

            ret
        } else {
            (0 .. nkpoints).collect::<Vec<usize>>()
        };

        let factor = val.factor.unwrap_or(1.0);

        let sel = Selection {
            label: label.to_string(),
            ispins,
            ikpoints,
            iatoms,
            iorbits,
            factor,
        };

        sel_vec.push(sel);
    }

    Ok(sel_vec)
}


#[derive(Clone, Serialize, Deserialize)]
pub enum SmearingMethod {
    Gaussian,
    Lorentz,
}


#[derive(Clone, Serialize, Deserialize)]
struct Configuration {
    #[serde(default = "Configuration::method_default")]
    method: SmearingMethod,

    #[serde(default = "Configuration::sigma_default")]
    sigma: f64,

    #[serde(default = "Configuration::procar_default")]
    procar: PathBuf,

    #[serde(default = "Configuration::outcar_default")]
    outcar: PathBuf,

    #[serde(default = "Configuration::txtout_default")]
    txtout: PathBuf,

    #[serde(default = "Configuration::htmlout_default")]
    htmlout: PathBuf,

    #[serde(default = "Configuration::notot_default")]
    notot: bool,

    pdos: Option<HashMap<String, RawSelection>>,
}

impl Configuration {
    pub fn method_default() -> SmearingMethod { SmearingMethod::Gaussian }
    pub fn sigma_default() -> f64 { 0.05 }
    pub fn procar_default() -> PathBuf { PathBuf::from("./PROCAR") }
    pub fn outcar_default() -> PathBuf { PathBuf::from("./OUTCAR") }
    pub fn txtout_default() -> PathBuf { PathBuf::from("./dos_raw.txt") }
    pub fn htmlout_default() -> PathBuf { PathBuf::from("./dos.html") }
    pub fn notot_default()  -> bool { false }
}



#[derive(Debug, StructOpt, Clone)]
#[structopt(setting = AppSettings::ColoredHelp,
            setting = AppSettings::ColorAuto)]
/// Calculate density of states from PROCAR and OUTCAR.
///
/// The fermi level is extracted from OUTCAR, and DOS is calculated by
/// smearing the band levels from PROCAR. The result may differ from DOSCAR.
pub struct Dos {
    #[structopt(short, long)]
    /// Projected DOS configuration file path.
    ///
    /// If left empty, only total DOS is calculated. The configuration template
    /// can be generated by `--gen-template` and then you can follow it.
    config: Option<PathBuf>,

    #[structopt(long)]
    /// Generate projected DOS configuration template.
    gen_template: bool,

    #[structopt(long, default_value = "./OUTCAR")]
    /// OUTCAR path
    outcar: PathBuf,

    #[structopt(long, default_value = "./PROCAR")]
    /// PROCAR path
    procar: PathBuf,

    #[structopt(long, default_value = "dos_raw.txt")]
    /// Save the raw data of projected DOS. Then you can replot it with more advanced tools.
    txtout: PathBuf,

    #[structopt(long, default_value = "dos.html")]
    /// Save the projected DOS plot as HTML. Then you can view it in the browser.
    ///
    /// Note: Your browser should be able to run plotly.js. Chrome, Safari, Edge, Firefox and
    /// etc. are supported.
    htmlout: PathBuf,

    #[structopt(long)]
    /// Print available orbits from PROCAR, this may be helpful when you write the configuration.
    show_orbits: bool
}


const TEMPLATE: &'static str = r#"# rsgrad DOS configuration in toml format.
method      = "Gaussian"        # smearing method
sigma       = 0.05              # smearing width, (eV)
procar      = "PROCAR"          # PROCAR path
outcar      = "OUTCAR"          # OUTCAR path
txtout      = "dos_raw.txt"     # save the raw data as "dos_raw.txt"
htmlout     = "dos.html"        # save the pdos plot as "dos.html"
notot       = false             # plot the total dos

[pdos.plot1]     # One label produces one plot, the labels CANNOT be repetitive.
spins   = "up down"         # for ISPIN = 2 system, "up" and "down" are available,
                            # for LSORBIT = .TRUE. system, "x" "y" "z" and "tot" are available.
kpoints = "1 3..7 -1"       # selects 1 3 4 5 6 7 and the last kpoint for pdos plot.
atoms   = "1 3..7 -1"       # selects 1 3 4 5 6 7 and the last atoms' projection for pdos plot.
orbits  = "s px dxy"         # selects the s px and dx orbits' projection for pdos plot.
factor  = 1.01               # the factor multiplied to this pdos

# The fields can be left blank, if you want select all the components for some fields,
# just comment them. You can comment fields with '#'
"#;


impl OptProcess for Dos {
    fn process(&self) -> Result<()> {
        if self.gen_template {
            let conf_filename = PathBuf::from("./dos.toml");

            info!("Generating selection dos configuration template ...");
            fs::write(&conf_filename, TEMPLATE.as_bytes())?;
            info!("Template file written to {:?}. Exiting", &conf_filename);
            
            return Ok(());
        }

        info!("Parsing PROCAR file {:?}", &self.procar);
        let procar  = Procar::from_file(&self.procar)?;
        let nlm     = procar.pdos.nlm.clone();
        let nkpts   = procar.kpoints.nkpoints as usize;

        if self.show_orbits {
            info!("Available orbits are\n {:?}", &nlm);
        }

        Ok(())
    }
}



#[cfg(test)]
mod test {
    use super::*;
    
    #[test]
    fn test_parse_rawconfig() {
        let nlm = "s     py     pz     px    dxy    dyz    dz2    dxz    dx2    tot" 
            .split_whitespace()
            .map(String::from)
            .collect::<Vec<_>>();
        let nions = 8usize;
        let nkpoints = 18usize;
        let nspin = 2;
        let is_ncl = false;

        let c: Configuration = toml::from_str(TEMPLATE).unwrap();
        let v = rawsel_to_sel(c.clone().pdos.unwrap(),
                              &nlm,
                              nions,
                              nkpoints,
                              nspin,
                              is_ncl).unwrap();
        assert_eq!(v[0].label, "plot1");
        assert_eq!(v[0].ispins, &[0, 1]);
        assert_eq!(v[0].iatoms, &[0, 2, 3, 4, 5, 6, 7]);
        assert_eq!(v[0].ikpoints, &[0, 2, 3, 4, 5, 6, 17]);
        assert_eq!(v[0].iorbits, &[0, 3, 4]);
        assert_eq!(v[0].factor, 1.01);

        let s = toml::to_string(&c).unwrap();
        println!("{}", s);
        println!("{:?}", v);
    }
}
