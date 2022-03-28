use std::{
    fs,
    path::PathBuf,
    time::Instant,
};

use indexmap::IndexMap;
use structopt::{
    StructOpt,
    clap::AppSettings,
};
use log::{
    info,
    debug,
};
use serde::{
    Serialize,
    Deserialize,
};
use anyhow::bail;
use toml;
use plotly;
use rayon::prelude::*;
use ndarray::{
    self,
    s,
};

use crate::{
    Result,
    OptProcess,
    Procar,
    //Outcar,
    vasp_parsers::outcar::GetEFermi,
    types::Vector,
    commands::common::{
        RawSelection,
        write_array_to_txt,
        CustomColor,
    }
};


const PI: f64 = std::f64::consts::PI;


#[derive(Clone, Debug)]
struct Selection {
    label:      String,
    ikpoints:   Vec<usize>,
    iatoms:     Vec<usize>,
    iorbits:    Vec<usize>,
    color:      Option<CustomColor>,
    factor:     f64,
}


fn rawsel_to_sel(r: IndexMap<String, RawSelection>, 
                 nlm: &[String], 
                 nions: usize, 
                 nkpoints: usize) -> Result<Vec<Selection>> {

    let mut sel_vec = vec![];

    for (label, val) in r.into_iter() {
        let ikpoints    = RawSelection::parse_ikpoints( val.kpoints.as_deref(), nkpoints)?;
        let iatoms      = RawSelection::parse_iatoms(   val.atoms.as_deref(),   nions)?;
        let iorbits     = RawSelection::parse_iorbits(  val.orbits.as_deref(),  nlm)?;
        let color       = if let Some(color) = val.color {
            Some( RawSelection::parse_color(&color)?)
        } else {
            None
        };
        let factor      = val.factor.unwrap_or(1.0);

        if factor < 0.0 { bail!("The factor cannot be negative."); }

        let sel = Selection {
            label: label.to_string(),
            ikpoints,
            iatoms,
            iorbits,
            color,
            factor,
        };

        sel_vec.push(sel);
    }

    Ok(sel_vec)
}


#[derive(Clone, Copy, Serialize, Deserialize, Debug)]
pub enum SmearingMethod {
    Gaussian,
    Lorentz,
}


#[derive(Clone, Serialize, Deserialize, Debug)]
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

    #[serde(default = "Configuration::totdos_default")]
    totdos: bool,

    #[serde(default = "Configuration::fill_default")]
    fill: bool,

    pdos: Option<IndexMap<String, RawSelection>>,
}

impl Configuration {
    pub fn method_default() -> SmearingMethod { SmearingMethod::Gaussian }
    pub fn sigma_default() -> f64 { 0.05 }
    pub fn procar_default() -> PathBuf { PathBuf::from("./PROCAR") }
    pub fn outcar_default() -> PathBuf { PathBuf::from("./OUTCAR") }
    pub fn txtout_default() -> PathBuf { PathBuf::from("./dos_raw.txt") }
    pub fn htmlout_default() -> PathBuf { PathBuf::from("./dos.html") }
    pub fn fill_default() -> bool { true }
    pub fn totdos_default()  -> bool { true }
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


impl Dos {
    // gaussian_smearing(x::AbstractArray, μ::Float64, σ=0.05) = @. exp(-(x-μ)^2 / (2*σ^2)) / (σ*sqrt(2π))
    fn smearing_gaussian(x: &[f64], mus: &[f64], sigma: f64, scales: &[f64]) -> Vector<f64> {
        let xlen = x.len();
        let clen = mus.len();
        let inv_two_sgm_sqr = 1.0 / (2.0 * sigma.powi(2));  // 1.0/(2*σ^2)
        let inv_sgm_sqrt2pi = 1.0 / (sigma * (2.0 * PI).sqrt()); // 1.0/(σ*sqrt(2π))

        let mut ret = Vector::<f64>::zeros(xlen);

        for c in 0 .. clen {
            ret.iter_mut()
                .zip(x.iter())
                .for_each(|(y, x)| {
                    *y += (-(x-mus[c]).powi(2) * inv_two_sgm_sqr).exp() * inv_sgm_sqrt2pi * scales[c];
                });
        }

        ret
    }

    // lorentz_smearing(x::AbstractArray, x0::Float64, Γ=0.05) = @. Γ/(2π) / ((x-x0)^2 + (Γ/2)^2)
    fn smearing_lorentz(x: &[f64], x0s: &[f64], gamma: f64, scales: &[f64]) -> Vector<f64> {
        let xlen = x.len();
        let clen = x0s.len();
        let gam_div_2pi = gamma / (2.0 * PI);  // Γ/(2π)
        let gam_half_sqr = (gamma / 2.0).powi(2); // (Γ/2)^2

        let mut ret = Vector::<f64>::zeros(xlen);

        for c in 0 .. clen {
            ret.iter_mut()
                .zip(x.iter())
                .for_each(|(y, x)| {
                    *y += gam_div_2pi / ((x - x0s[c]).powi(2) + gam_half_sqr) * scales[c];
                })
        }

        ret
    }

    fn apply_smearing(x: &[f64], centers: &[f64], width: f64, method: SmearingMethod, scales: Option<&[f64]>) -> Vector<f64> {
        let clen = centers.len();
        let mut fac = vec![1.0; 0];

        let scales = if let Some(factors) = scales {
            assert_eq!(factors.len(), clen, "[DOSsmear]: factors' length inconsistent with smearing peaks");
            factors
        } else {
            fac.resize(clen, 1.0);
            &fac
        };

        match method {
            SmearingMethod::Lorentz  => Self::smearing_lorentz(x, centers, width, scales),
            SmearingMethod::Gaussian => Self::smearing_gaussian(x, centers, width, scales),
        }
    }


    fn gen_totdos(xvals: &[f64], procar: &Procar, sigma: f64, method: SmearingMethod) -> Vector<f64> {
        let nspin       = procar.pdos.nspin as usize;
        let nkpoints    = procar.pdos.nkpoints as usize;
        let mut totdos  = vec![];

        let norm = procar.kpoints.weights.sum();
        let weights = &procar.kpoints.weights / norm;

        for ispin in 0 .. nspin {
            let mut tdos = Vector::<f64>::zeros(xvals.len());
            for ikpoint in 0 .. nkpoints {
                let eigs = procar.pdos.eigvals.slice(s![ispin, ikpoint, ..]).to_slice().unwrap();
                if 0 == ispin {
                    tdos += &(Self::apply_smearing(xvals, eigs, sigma, method, None) * weights[ikpoint]);
                } else {
                    tdos -= &(Self::apply_smearing(xvals, eigs, sigma, method, None) * weights[ikpoint]);
                }
            }

            let tdos = if 0 == ispin {
                tdos.into_raw_vec()
            } else {
                let mut r = tdos.into_raw_vec();
                r.reverse();
                r
            };

            totdos.push(tdos);
        }

        totdos.into_iter()
            .map(|x| x.into_iter())
            .flatten()
            .collect()
    }

    fn gen_pdos(xvals: &[f64], procar: &Procar, selection: &Selection, sigma: f64, method: SmearingMethod) -> Vector<f64> {
        let nspin       = procar.pdos.nspin as usize;
        let nbands      = procar.pdos.nbands as usize;
        let factor      = selection.factor;
        let mut totdos  = Vec::<Vec<f64>>::new();

        let norm        = procar.kpoints.weights.sum();
        let kptweights  = &procar.kpoints.weights / norm;

        // pdos.projected = [ispin, ikpoint, iband, iion, iorbit]

        for ispin in 0 .. nspin {
            let mut tdos = Vector::<f64>::zeros(xvals.len());
            for ikpoint in selection.ikpoints.iter().copied() {
                let eigs = procar.pdos.eigvals.slice(s![ispin, ikpoint, ..]).to_slice().unwrap();
                let bandweights = (0 .. nbands)
                    .into_iter()
                    .map(|iband| {
                        let mut wht = 0.0;
                        for iion in selection.iatoms.iter().copied() {
                            for iorbit in selection.iorbits.iter().copied() {
                                wht += procar.pdos.projected[[ispin, ikpoint, iband, iion, iorbit]];
                            }
                        }
                        wht
                    }).collect::<Vec<f64>>();

                if 0 == ispin {
                    tdos += &(Self::apply_smearing(xvals, eigs, sigma, method, Some(&bandweights)) * kptweights[ikpoint]);
                } else {
                    tdos -= &(Self::apply_smearing(xvals, eigs, sigma, method, Some(&bandweights)) * kptweights[ikpoint]);
                }

            }

            tdos *= factor;

            let tdos = if 0 == ispin {
                tdos.into_raw_vec()
            } else {
                let mut r = tdos.into_raw_vec();
                r.reverse();
                r
            };

            totdos.push(tdos);
        }

        totdos.into_iter()
            .map(|x| x.into_iter())
            .flatten()
            .collect()
    }
}


const TEMPLATE: &'static str = r#"# rsgrad DOS configuration in toml format.
# multiple tokens inside string are seperated by whitespace
method      = "Gaussian"        # smearing method
sigma       = 0.05              # smearing width, (eV)
procar      = "PROCAR"          # PROCAR path
outcar      = "OUTCAR"          # OUTCAR path
txtout      = "dos_raw.txt"     # save the raw data as "dos_raw.txt"
htmlout     = "dos.html"        # save the pdos plot as "dos.html"
totdos      = true              # plot the total dos
fill        = true              # fill the plot to x axis or not

[pdos.plot1]                  # One label produces one plot, the labels CANNOT be repetitive.
                              # each the label is 'plot1', to add more pdos, write '[pdos.plot2]' and so on.
kpoints = "1 3..7 -1"         # selects 1 3 4 5 6 7 and the last kpoint for pdos plot.
atoms   = "1 3..7 -1"         # selects 1 3 4 5 6 7 and the last atoms' projection for pdos plot.
orbits  = "s px dxy"          # selects the s px and dx orbits' projection for pdos plot.
factor  = 1.01                # the factor multiplied to this pdos

[pdos.plot2]                  # One label produces one plot, the labels CANNOT be repetitive.
                              # each the label is 'plot1', to add more pdos, write '[pdos.plot2]' and so on.
kpoints = "1 3..7 -1"         # selects 1 3 4 5 6 7 and the last kpoint for pdos plot.
atoms   = "1 3..7 -1"         # selects 1 3 4 5 6 7 and the last atoms' projection for pdos plot.
orbits  = "s px dxy"          # selects the s px and dx orbits' projection for pdos plot.
factor  = 1.01                # the factor multiplied to this pdos

# The fields can be left blank, if you want select all the components for some fields,
# just comment them. You can comment fields with '#'
"#;


impl OptProcess for Dos {
    fn process(&self) -> Result<()> {
        if self.gen_template {
            let conf_filename = PathBuf::from("./pdos.toml");

            info!("Generating selection dos configuration template ...");
            fs::write(&conf_filename, TEMPLATE.as_bytes())?;
            info!("Template file written to {:?}. Exiting", &conf_filename);
            
            return Ok(());
        }

        let config = if let Some(config) = self.config.as_ref() {
            info!("Reading PDOS configuration from {:?}", self.config.as_ref());
            let config = fs::read_to_string(config)?;
            let config: Configuration = toml::from_str(&config)?;

            debug!("{:#?}", &config);

            Some(config)
        } else {
            None
        };

        let procar     = if let Some(cfg) = config.as_ref() {  &cfg.procar } else {  &self.procar };
        let outcar     = if let Some(cfg) = config.as_ref() {  &cfg.outcar } else {  &self.outcar };
        let txtout     = if let Some(cfg) = config.as_ref() {  &cfg.txtout } else {  &self.txtout };
        let htmlout    = if let Some(cfg) = config.as_ref() { &cfg.htmlout } else { &self.htmlout };
        let sigma      = if let Some(cfg) = config.as_ref() {    cfg.sigma } else {          0.05 };
        let method     = if let Some(cfg) = config.as_ref() {   cfg.method } else { SmearingMethod::Gaussian };
        let is_totdos  = if let Some(cfg) = config.as_ref() {   cfg.totdos } else {          true };
        let to_fill    = if let Some(cfg) = config.as_ref() {     cfg.fill } else {          true };

        if sigma < 0.0 {
            bail!("[DOS]: sigma cannot be negative.");
        }

        info!("Parsing PROCAR file {:?}", procar);
        let mut procar  = Procar::from_file(procar)?;
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
        info!("Found Fermi level = {}, eigenvalues will be shifted.", efermi);

        let selections = if config.as_ref().is_some() {
            if let Some(pdos) = config.clone().unwrap().pdos {
                Some(rawsel_to_sel(pdos, &nlm, nions, nkpts)?)
            } else {
                None
            }
        } else {
            None
        };

        procar.pdos.eigvals -= efermi;

        let emin = procar.pdos.eigvals
            .iter()
            .cloned()
            .reduce(f64::min)
            .unwrap();
        let emax = procar.pdos.eigvals
            .iter()
            .cloned()
            .reduce(f64::max)
            .unwrap();


        let mut plot = plotly::Plot::new();

        let nedos = (emax - emin).ceil() as usize * 400;  // 400 points per eV
        let xvals = Vector::<f64>::linspace(emin-2.0, emax+2.0, nedos);
        let xvals_plot = if nspin == 1 {
            xvals.clone()
        } else {
            xvals.clone()
                .into_iter()
                .chain(xvals.clone().into_raw_vec().into_iter().rev())
                .collect::<Vector<f64>>()
        };

        let totdos = Self::gen_totdos(xvals.as_slice().unwrap(), &procar, sigma, method);

        let mut labels = vec!["E-Ef", "TotDOS"]
            .into_iter()
            .map(String::from)
            .collect::<Vec<_>>();

        let mut raw_dats = vec![
            xvals_plot.clone(), 
            totdos.clone()
        ];

        let fill_type = if to_fill {
            plotly::common::Fill::ToZeroY
        } else {
            plotly::common::Fill::None
        };

        if is_totdos {
            info!("Plotting Total DOS ...");
            let now = Instant::now();
            let tr = plotly::Scatter::from_array(xvals_plot.clone(), totdos.clone())
                .mode(plotly::common::Mode::Lines)
                .marker(plotly::common::Marker::new()
                        .color(plotly::NamedColor::Black))
                .fill(fill_type.clone())
                .name("Total DOS");
            plot.add_trace(tr);
            info!("Total DOS plot time usage: {:?}", now.elapsed());
        };

        if let Some(selections) = selections {
            info!("Plotting PDOS ...");
            let now = Instant::now();

            let doses = selections
                .into_par_iter()
                .map(|sel| {
                    (
                        Self::gen_pdos(xvals.as_slice().unwrap(), &procar, &sel, sigma, method),
                        sel.label,
                        sel.color,
                    )
                })
                .collect::<Vec<_>>();

            for (dos, label, color) in doses.into_iter() {
                let mut marker = plotly::common::Marker::new();
                if let Some(c) = color {
                    marker = marker.color(c);
                };

                let tr = plotly::Scatter::from_array(xvals_plot.clone(), dos.clone())
                    .mode(plotly::common::Mode::Lines)
                    .marker(marker)
                    .fill(fill_type.clone())
                    .name(&label);
                labels.push(label);
                raw_dats.push(dos);
                plot.add_trace(tr);
            }

            info!("PDOS plot time usage: {:?}", now.elapsed());
        }


        plot.use_local_plotly();
        let layout = plotly::Layout::new()
            .title(plotly::common::Title::new("Density of States"))
            .y_axis(plotly::layout::Axis::new()
                    .title(plotly::common::Title::new("DOS (arb. unit)"))
                    .zero_line(true)
                    .fixed_range(false)
                    )
            .x_axis(plotly::layout::Axis::new()
                    .title(plotly::common::Title::new("E-Ef (eV)"))
                    .range(vec![-1.0, 6.0])
                    .zero_line(true)
                    .range_slider(plotly::layout::RangeSlider::new().visible(true))
                    );
        plot.set_layout(layout);

        info!("Writing DOS plot to {:?}", htmlout);
        plot.to_html(&htmlout);

        info!("Writing raw DOS data to {:?}", txtout);
        let label = labels.join(" ");
        let raw_dats = raw_dats.iter().map(|x| x).collect::<Vec<_>>();
        write_array_to_txt(txtout, raw_dats, &label)?;

        Ok(())
    }
}



#[cfg(test)]
mod test {
    use super::*;
    
    #[test]
    fn test_parse_rawconfig() {
        let nlm = "s     py     pz     px    dxy    dyz    dz2    dxz    dx2" 
            .split_whitespace()
            .map(String::from)
            .collect::<Vec<_>>();
        let nions = 8usize;
        let nkpoints = 18usize;

        let c: Configuration = toml::from_str(TEMPLATE).unwrap();
        let v = rawsel_to_sel(c.clone().pdos.unwrap(),
                              &nlm,
                              nions,
                              nkpoints).unwrap();
        assert_eq!(v[0].label, "plot1");
        assert_eq!(v[0].iatoms, &[0, 2, 3, 4, 5, 6, 7]);
        assert_eq!(v[0].ikpoints, &[0, 2, 3, 4, 5, 6, 17]);
        assert_eq!(v[0].iorbits, &[0, 3, 4]);
        assert_eq!(v[0].factor, 1.01);

        let s = toml::to_string(&c).unwrap();
        println!("{}", s);
        println!("{:?}", v);
    }

}
