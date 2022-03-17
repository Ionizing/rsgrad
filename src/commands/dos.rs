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
};
use anyhow::{
    //anyhow,
    Context,
    //Error,
    bail,
};
use toml;

use crate::{
    Result,
    OptProcess,
    Procar,
    //Outcar,
    vasp_parsers::outcar::GetEFermi,
    types::{
        range_parse,
        index_transform,
        Vector,
    },
};


const PI: f64 = std::f64::consts::PI;


#[derive(Clone, Serialize, Deserialize)]
struct RawSelection {
    spins:      Option<String>,
    kpoints:    Option<String>,
    atoms:      Option<String>,
    orbits:     Option<String>,
    factor:     Option<f64>,
}


impl RawSelection {
    fn parse_ispins(input: Option<&str>, nspin: usize, is_ncl: bool) -> Result<Vec<usize>> {
        if let Some(spins) = input {
            if is_ncl {
                spins.split_whitespace()
                    .map(|x| match x.to_ascii_lowercase().as_ref() {
                        "x"         => Ok(1usize),
                        "y"         => Ok(2usize),
                        "z"         => Ok(3usize),
                        "t" | "tot" => Ok(0usize),
                        _ =>
bail!("[DOS]: Invalid spin component selected: `{}`, available components are `x`, `y`, `z` and `tot`", x)
                    }).collect::<Result<Vec<_>>>()
            } else if nspin == 2 {
                spins.split_whitespace()
                    .map(|x| match x.to_ascii_lowercase().as_ref() {
                        "u" | "up"           => Ok(0usize),
                        "d" | "dn" | "down"  => Ok(1usize),
                        _ =>
bail!("[DOS]: Invalid spin component selected: `{}`, available components are `up` and `down`", x)
                    }).collect::<Result<Vec<_>>>()
            } else {
                warn!("[DOS]: This system is not spin-polarized, only one spin component is available, selected by default.");
                Ok(vec![0usize])
            }
        } else {
            Ok(if is_ncl {
                vec![0usize]  // 'tot' part
            } else if nspin == 2 {
                vec![0usize, 1usize]  // spin up and spin down
            } else {
                vec![0usize]  // only one spin available
            })
        }
    }

    fn parse_iatoms(input: Option<&str>, nions: usize) -> Result<Vec<usize>> {
        if let Some(atoms) = input {
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
            Ok(ret)
        } else {  // All atoms selected if left blank
            Ok((0 .. nions).collect::<Vec<usize>>())
        }
    }

    fn parse_ikpoints(input: Option<&str>, nkpoints: usize) -> Result<Vec<usize>> {
        if let Some(kpoints) = input {
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
            Ok(ret)
        } else {
            Ok((0 .. nkpoints).collect::<Vec<usize>>())
        }
    }

    fn parse_iorbits(input: Option<&str>, nlm: &[String]) -> Result<Vec<usize>> {
        if let Some(orbits) = input {
            orbits.split_whitespace()
                .map(|x| {
                    nlm.iter().position(|x2| x2 == &x)
                        .context(format!("Selected orbit {:?} not available in {:?}", x, &nlm))
                })
            .collect::<Result<Vec<_>>>()
        } else {
            Ok((0 .. nlm.len()).collect::<Vec<usize>>())
        }
    }
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
        let ispins = RawSelection::parse_ispins(    val.spins.as_deref(),   nspin,  is_ncl)?;
        let ikpoints = RawSelection::parse_ikpoints(val.kpoints.as_deref(), nkpoints)?;
        let iatoms = RawSelection::parse_iatoms(    val.atoms.as_deref(),   nions)?;
        let iorbits = RawSelection::parse_iorbits(  val.orbits.as_deref(),  nlm)?;
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


#[derive(Clone, Copy, Serialize, Deserialize)]
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
    /// Print brief info of PROCAR, this may be helpful when you write the configuration.
    show_brief: bool
}


// gaussian_smearing(x::AbstractArray, μ::Float64, σ=0.05) = @. exp(-(x-μ)^2 / (2*σ^2)) / (σ*sqrt(2π))
fn smearing_gaussian(x: &[f64], mus: &[f64], sigma: f64) -> Vector<f64> {
    let len = x.len();
    let inv_two_sgm_sqr = 1.0 / (2.0 * sigma.powi(2));  // 1.0/(2*σ^2)
    let inv_sgm_sqrt2pi = 1.0 / (sigma * (2.0 * PI).sqrt()); // 1.0/(σ*sqrt(2π))

    let mut ret = Vector::<f64>::zeros(len);

    for i in 0..len {
        ret[i] += mus.as_ref().iter()
            .map(|c| (-(x[i]-c).powi(2) * inv_two_sgm_sqr).exp() * inv_sgm_sqrt2pi)
            .sum::<f64>();
    }

    ret
}

// lorentz_smearing(x::AbstractArray, x0::Float64, Γ=0.05) = @. Γ/(2π) / ((x-x0)^2 + (Γ/2)^2)
fn smearing_lorentz(x: &[f64], x0s: &[f64], gamma: f64) -> Vector<f64> {
    let len = x.len();
    let gam_div_2pi = gamma / (2.0 * PI);  // Γ/(2π)
    let gam_half_sqr = (gamma / 2.0).powi(2); // (Γ/2)^2

    let mut ret = Vector::<f64>::zeros(len);

    for i in 0..len {
        ret[i] += x0s.as_ref().iter()
            .map(|c| gam_div_2pi / ((x[i] - c).powi(2) + gam_half_sqr))
            .sum::<f64>();
    }

    ret
}

fn apply_smearing(x: &[f64], centers: &[f64], width: f64, method: SmearingMethod) -> Vector<f64> {
    match method {
        SmearingMethod::Lorentz  => smearing_lorentz(x, centers, width),
        SmearingMethod::Gaussian => smearing_gaussian(x, centers, width),
    }
}


fn gen_totdos() {
    todo!()
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
kpoints = "1 3..7 -1"       # selects 1 3 4 5 6 7 and the last kpoint for pdos plot, starts from 1.
atoms   = "1 3..7 -1"       # selects 1 3 4 5 6 7 and the last atoms' projection for pdos plot, starts from 1.
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

        let config = if let Some(config) = self.config.as_ref() {
            info!("Reading PDOS configuration from {:?}", self.config.as_ref());
            let config = fs::read_to_string(config)?;
            let config: Configuration = toml::from_str(&config)?;
            Some(config)
        } else {
            None
        };

        let procar  = if let Some(cfg) = config.as_ref() {  &cfg.procar } else {  &self.procar };
        let outcar  = if let Some(cfg) = config.as_ref() {  &cfg.outcar } else {  &self.outcar };
        let txtout  = if let Some(cfg) = config.as_ref() {  &cfg.txtout } else {  &self.txtout };
        let htmlout = if let Some(cfg) = config.as_ref() { &cfg.htmlout } else { &self.htmlout };
        let sigma   = if let Some(cfg) = config.as_ref() {    cfg.sigma } else {          0.05 };
        let method  = if let Some(cfg) = config.as_ref() {   cfg.method } else { SmearingMethod::Gaussian };
        let notot   = if let Some(cfg) = config.as_ref() {    cfg.notot } else {         false };
        

        info!("Parsing PROCAR file {:?}", procar);
        let procar  = Procar::from_file(procar)?;
        let nlm     = procar.pdos.nlm.clone();
        let nkpts   = procar.kpoints.nkpoints as usize;
        let nions   = procar.pdos.nions as usize;
        let is_ncl  = procar.pdos.lsorbit;
        let nspin   = procar.pdos.nspin as usize;
        let nbands  = procar.pdos.nbands as usize;

        if self.show_brief {
            info!("Brief info of current PROCAR:
  orbitals = {:?}
  NKPTS  = {}
  NIONS  = {}
  IS_NCL = {}
  NSPIN  = {}
  NBANDS = {}", &nlm, nkpts, nions, is_ncl, nspin, nbands);
            return Ok(());
        }

        info!("Parsing OUTCAR file {:?} for Fermi level", outcar);
        let efermi = fs::read_to_string(outcar)?.get_efermi()?;
        info!("Found Fermi level = {}", efermi);

        Ok(())
    }
}



#[cfg(test)]
mod test {
    use super::*;

    const TEST_TEMPLATE: &'static str = r#"# rsgrad DOS configuration in toml format.
method      = "Gaussian"        # smearing method
sigma       = 0.05              # smearing width, (eV)
procar      = "PROCAR"          # PROCAR path
outcar      = "OUTCAR"          # OUTCAR path
txtout      = "dos_raw.txt"     # save the raw data as "dos_raw.txt"
htmlout     = "dos.html"        # save the pdos plot as "dos.html"
notot       = false             # plot the total dos

[pdos.plot1]                  # One label produces one plot, the labels CANNOT be repetitive.
                              # each the label is 'plot1', to add more pdos, write '[pdos.plot2]' and so on.
spins   = "up down"           # for ISPIN = 2 system, "up" and "down" are available,
                              # for LSORBIT = .TRUE. system, "x" "y" "z" and "tot" are available.
kpoints = "1 3..7 -1"         # selects 1 3 4 5 6 7 and the last kpoint for pdos plot.
atoms   = "1 3..7 -1"         # selects 1 3 4 5 6 7 and the last atoms' projection for pdos plot.
orbits  = "s px dxy"          # selects the s px and dx orbits' projection for pdos plot.
factor  = 1.01                # the factor multiplied to this pdos

# The fields can be left blank, if you want select all the components for some fields,
# just comment them. You can comment fields with '#'
"#;
    
    #[test]
    fn test_parse_rawconfig() {
        let nlm = "s     py     pz     px    dxy    dyz    dz2    dxz    dx2" 
            .split_whitespace()
            .map(String::from)
            .collect::<Vec<_>>();
        let nions = 8usize;
        let nkpoints = 18usize;
        let nspin = 2;
        let is_ncl = false;

        let c: Configuration = toml::from_str(TEST_TEMPLATE).unwrap();
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

    #[test]
    fn test_parse_rawconfig_2() {
        let nions = 8usize;
        let nkpoints = 18usize;
        let nspin = 1;
        let is_ncl = true;
        let nlm = "s     py     pz     px    dxy    dyz    dz2    dxz    dx2" 
            .split_whitespace()
            .map(String::from)
            .collect::<Vec<_>>();
        let cfg = r#"[pdos.plot1]
spins = "down"
"#;
        let c:Configuration = toml::from_str(cfg).unwrap();
        let v = rawsel_to_sel(c.clone().pdos.unwrap(),
                              &nlm,
                              nions,
                              nkpoints,
                              nspin,
                              is_ncl);
        assert!(v.is_err());

        let nions = 8usize;
        let nkpoints = 18usize;
        let nspin = 1;
        let is_ncl = false;
        let nlm = "s     py     pz     px    dxy    dyz    dz2    dxz    dx2" 
            .split_whitespace()
            .map(String::from)
            .collect::<Vec<_>>();
        let cfg = r#"[pdos.plot1]
spins = "down"
"#;
        
    }
}
