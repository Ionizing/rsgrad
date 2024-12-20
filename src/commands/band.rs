use std::{
    fs,
    path::PathBuf,
    time::Instant,
};

use indexmap::IndexMap;
use clap::Args;
use log::{
    info,
    warn,
    debug,
};
use rayon::{
    self,
    prelude::*,
};
use anyhow::{
    bail,
    anyhow,
    Context,
};
use serde::{
    Serialize,
    Deserialize,
    Deserializer,
};
use ndarray::{
    s,
    arr2,
    Array5,
    concatenate,
};
use plotly::{
    self,
    common::color::NamedColor,
    Plot,
    Scatter,
    common::ColorScalePalette,
    layout::{
        Shape,
        ShapeType,
        ShapeLine,
        Layout,
        ItemSizing,
    }
};

use crate::{
    Result,
    OptProcess,
    Procar,
    Outcar,
    Poscar,
    types::{
        Vector,
        Matrix,
        Cube,
        Axis,
    },
    commands::common::{
        write_array_to_txt,
        RawSelection,
        generate_plotly_configuration,
    }
};


const THRESHOLD: f64 = 1E-6;


#[derive(Clone, Debug)]
struct Selection {
    label:      String,
    ispins:     Vec<usize>,
    iatoms:     Vec<usize>,
    iorbits:    Vec<usize>,
    color:      Option<String>,
}


fn rawsel_to_sel(r: IndexMap<String, RawSelection>, 
                 nspin: usize,
                 is_ncl: bool,
                 nlm: &[String], 
                 nions: usize) -> Result<Vec<Selection>> {

    let mut sel_vec = vec![];

    for (label, val) in r.into_iter() {
        let ispins      = RawSelection::parse_ispins(   val.spins.as_deref(),   nspin, is_ncl)?;
        let iatoms      = RawSelection::parse_iatoms(   val.atoms.as_deref(),   nions)?;
        let iorbits     = RawSelection::parse_iorbits(  val.orbits.as_deref(),  nlm)?;
        let color       = if let Some(color) = val.color {
            Some(RawSelection::parse_color(&color)?)
        } else {
            None
        };

        let sel = Selection {
            label: label.to_string(),
            ispins,
            iatoms,
            iorbits,
            color,
        };

        sel_vec.push(sel);
    }

    Ok(sel_vec)
}


#[derive(Clone, Serialize, Deserialize, Debug)]
#[serde(rename_all = "kebab-case",
        deny_unknown_fields)]
struct Configuration {
    kpoint_labels: Option<Vec<String>>,

    #[serde(default = "Configuration::procar_default")]
    procar: PathBuf,

    #[serde(default = "Configuration::outcar_default")]
    outcar: PathBuf,

    #[serde(default = "Configuration::txtout_prefix_default")]
    txtout_prefix: String,

    #[serde(default = "Configuration::htmlout_default")]
    htmlout: PathBuf,

    segment_ranges: Option<Vec<(usize, usize)>>,

    ncl_spinor: Option<Axis>,

    #[serde(default = "Configuration::colormap_default",
            deserialize_with = "Configuration::colormap_de")]
    colormap: ColorScalePalette,

    efermi: Option<f64>,

    #[serde(default = "Configuration::ylim_default")]
    ylim: (f64, f64),

    pband: Option<IndexMap<String, RawSelection>>,
}

impl Configuration {
    pub fn procar_default()         -> PathBuf { PathBuf::from("./PROCAR") }
    pub fn outcar_default()         -> PathBuf { PathBuf::from("./OUTCAR") }
    pub fn txtout_prefix_default()  -> String  { String::from("./band_raw") }
    pub fn htmlout_default()        -> PathBuf { PathBuf::from("./band.html") }
    pub fn colormap_default()       -> ColorScalePalette { ColorScalePalette::Jet }
    pub fn colormap_de<'de, D: Deserializer<'de>>(d: D) -> std::result::Result<ColorScalePalette, D::Error> {
        let s: Option<String> = Deserialize::deserialize(d)?;
        if let Some(s) = s {
            let cmap = RawSelection::parse_colormap(&s);
            match cmap {
                Ok(c) => { Ok(c) },
                Err(e) => { Err(serde::de::Error::custom(e.to_string())) },
            }
        } else {
            Ok(ColorScalePalette::Jet)
        }
    }
    pub fn ylim_default()            -> (f64, f64) { (-1.0, 6.0) }
}


#[derive(Debug, Clone, Args)]
#[command(allow_negative_numbers = true)]
/// Plot bandstructure and projected bandstructure
/// 
/// You may need to generate a template and then use the `--config` option to plot the projected
/// bandstructure.
pub struct Band {
    #[arg(short, long)]
    /// Band structure plot configuration file path.
    ///
    /// If left empty, only bare band is calculated. The configuration template
    /// can be generated by `--gen-template` and then you can follow it.
    config: Option<PathBuf>,

    #[arg(long)]
    /// Generate band structure plot configuration template.
    gen_template: bool,

    #[arg(long)]
    /// Set the E-fermi given from SCF's OUTCAR and it will be the reference energy level.
    ///
    /// If it is unset, the fermi level would be read from OUTCAR of non-scf calculation,
    /// i.e. bandstructure calculation, which may be a little different from scf's.
    efermi: Option<f64>,

    #[arg(short, long, num_args(0..))]
    /// Symbols for high symmetry points on the kpoint path.
    kpoint_labels: Option<Vec<String>>,

    #[arg(long, default_value = "./PROCAR")]
    /// PROCAR path.
    ///
    /// The band level and projected band info are extracted from PROCAR.
    procar: PathBuf,

    #[arg(long, default_value = "./OUTCAR")]
    /// OUTCAR path.
    ///
    /// The fermi level and lattice info are extracted from OUTCAR.
    outcar: PathBuf,

    #[arg(long, default_value = "jet", value_parser(RawSelection::parse_colormap))]
    colormap: ColorScalePalette,

    #[arg(long, value_enum, ignore_case = true)]
    ncl_spinor: Option<Axis>,

    #[arg(long, default_value = "band_raw")]
    /// Save the raw data of band structure.
    ///
    /// Then you can replot it with more advanced tools.
    txtout_prefix: String,

    #[arg(long, default_value = "band.html")]
    /// Save the band structure plot as HTML.
    ///
    /// Note: Your browser should be able to run plotly.js. Chrome, Safari, Edge, Firefox and
    /// etc. are supported.
    htmlout: PathBuf,

    #[arg(long)]
    /// Open the browser and show the plot immediately.
    show: bool,

    #[arg(long)]
    /// Render the plot and print the rendered code to stdout.
    to_inline_html: bool,

    #[arg(long, default_values = &["-1", "6"], num_args(2))]
    /// Set the y-range of the plot.
    ylim: Vec<f64>,
}


// Extra TODO: shift eigvals to E-fermi
impl Band {
    /// Given a `nkpoints*3` matrix to generate k-point path for bandstructure
    /// `segment_ranges` are the start and end indices of each, closed interval.
    /// the range can be in reversed direction. Index starts from 1.
    fn gen_kpath(kpoints: &Matrix<f64>, bcell: &[[f64; 3]; 3], segment_ranges: &[(usize, usize)])-> (Vector<f64>, Vec<f64>) {
        let bcell = arr2(bcell);
        let kpoints = kpoints.dot(&bcell);

        let kpath = segment_ranges
            .iter()
            .cloned()
            .flat_map(|(beg, end)| {
                let (rng, lrng) = if beg < end {    // Forward range (e.g. 1..5)
                    (s![beg..end,    ..], s![(beg-1)..(end-1),    ..])
                } else {                            // Backward range (e.g.  5..1)
                    (s![end..beg;-1, ..], s![(end-1)..(beg-1);-1, ..])
                };

                let mut kdiffs = (kpoints.slice(rng).to_owned() - kpoints.slice(lrng))
                    .outer_iter()
                    .map(|v| f64::sqrt(v[0]*v[0] + v[1]*v[1] + v[2]*v[2]))
                    .collect::<Vec<f64>>();

                kdiffs.insert(0, 0.0);
                kdiffs.into_iter()
            })
            .scan(0.0, |acc, x| {  // equal to cumsum
                *acc += x;
                Some(*acc)
            })
            .collect::<Vector<f64>>();

        let mut kxs = segment_ranges.iter()
            .map(|(ki, kj)| {  // ki and kj form a closed interval
                if ki > kj {
                    ki - kj + 1
                } else {
                    kj - ki + 1
                }
            })
            .scan(0, |acc, x| {     // equal to cumsum
                *acc += x;
                Some(*acc - 1)
            })
            .map(|ik| kpath[ik])
            .collect::<Vec<f64>>();

        kxs.insert(0, 0.0);

        (kpath, kxs)
    }

    /// Try to partition the path according to the equality of k-points
    fn find_segments(kpoints: &Matrix<f64>) -> Result<Vec<(usize, usize)>> {
        let len = kpoints.shape()[0];

        if len <= 3 {
            bail!("Too few k-points, consider add more k-points in calculation.");
        }

        let mut boundaries = Vec::<(usize, usize)>::new();
        let mut last_bound = 1usize;

        for i in 1 .. len {
            if (kpoints.index_axis(ndarray::Axis(0), i).to_owned() - 
                kpoints.index_axis(ndarray::Axis(0), i-1))
                .iter().all(|d| d.abs() < THRESHOLD) {  // if kpt[i, :] == kpt[i-1, :]

                boundaries.push((last_bound, i));
                last_bound = i + 1;
            }
        }

        boundaries.push((last_bound, len));
        Ok(boundaries)
    }


    /// Return: [ispin, ikpoint, iband]
    fn gen_rawband(eigvals: &Cube<f64>, segment_ranges: &[(usize, usize)]) -> Cube<f64> {
        let bands = segment_ranges.iter()
            .cloned()
            .map(|(beg, end)| {
                let rng = if beg < end {
                    s![.., beg-1 .. end, ..]
                } else {
                    s![.., end-1 .. beg;-1, ..]
                };
                eigvals.slice(rng)
            })
            .collect::<Vec<_>>();
        concatenate(ndarray::Axis(1), &bands).unwrap()
    }

    /// Return: [spins, kpoints, bands, nions, nnlm]
    fn gen_cropped_projections(proj: &Array5<f64>, segment_ranges: &[(usize, usize)]) -> Array5<f64> {
        let projections = segment_ranges.iter()
            .cloned()
            .map(|(beg, end)| {
                let rng = if beg < end {
                    s![.., beg-1 ..end, .., .., ..]
                } else {
                    s![.., end-1 ..beg;-1, .., .., ..]
                };
                proj.slice(rng)
            })
            .collect::<Vec<_>>();
        concatenate(ndarray::Axis(1), &projections).unwrap()
    }

    /// Plot the band dispersion only
    fn plot_rawband(plot: &mut Plot, kpath: Vector<f64>, cropped_eigvals: &Cube<f64>) {
        let nspin     = cropped_eigvals.shape()[0];
        let nkpoints  = cropped_eigvals.shape()[1];
        let nbands    = cropped_eigvals.shape()[2];

        assert_eq!(kpath.len(), nkpoints);    // cropped_eigvals[ispin, ikpoint, iband]

        let getcolor = |ispin: usize| {
            match (nspin, ispin) {
                (1, _) => NamedColor::Black,
                (2, 0) => NamedColor::Red,
                (2, 1) => NamedColor::Blue,
                _ => unreachable!("Invalid spin index"),
            }
        };

        for ispin in 0 .. nspin {
            (0 .. nbands)
                .for_each(|iband| {
                    let dispersion = cropped_eigvals.slice(s![ispin, .., iband]).to_owned();
                    let show_legend = 0 == iband;
                    let legend_name = match (nspin, ispin) {
                        (1, _) => "Band Dispersion",
                        (2, 0) => "Spin Up",
                        (2, 1) => "Spin Down",
                        _ => unreachable!("Invalied spin index"),
                    };

                    let hover_template0 = match (nspin, ispin) {
                        (1, _) => format!("Band#: {}<br>", iband + 1),
                        (2, 0) => format!("Spin up<br>Band#: {}<br>", iband + 1),
                        (2, 1) => format!("Spin Down<br>Band#: {}<br>", iband + 1),
                        _ => unreachable!("Only two spin components available"),
                    };

                    let tr = Scatter::from_array(kpath.clone(), dispersion)
                        .mode(plotly::common::Mode::Lines)
                        .marker(plotly::common::Marker::new().color(getcolor(ispin)))
                        .legend_group("Total bandstructure")
                        .show_legend(show_legend)
                        .hover_info(plotly::common::HoverInfo::Text)
                        .hover_template(hover_template0 + "E-Ef: %{y:.4f} eV")
                        .name(legend_name);
                    plot.add_trace(tr);
                });
        }
    }

    fn plot_boundaries(layout: &mut Layout, kxs: &[f64]) {
        kxs.iter()
            .cloned()
            .for_each(|k| {         // add vlines to canvas to identify high-symmetry points
                let shape = Shape::new()
                    .shape_type(ShapeType::Line)
                    .x0(k).y0(0.0)
                    .x1(k).y1(1.0)
                    .x_ref("x").y_ref("paper")
                    .line(ShapeLine::new()
                          .color(NamedColor::Black)
                          .width(0.7));
                layout.add_shape(shape);
            });

        let kmax = kxs.iter().last().cloned().unwrap();

        layout.add_shape(
            Shape::new()
            .shape_type(ShapeType::Line)
            .x0(0.0).y0(0.0)
            .x1(kmax).y1(0.0)       // add hline at the bottom
            .x_ref("x")
            .y_ref("paper")
            .line(ShapeLine::new()
                  .color(NamedColor::Black)
                  .width(0.7))
            );
        layout.add_shape(
            Shape::new()
            .shape_type(ShapeType::Line)
            .x0(0.0).y0(1.0)
            .x1(kmax).y1(1.0)       // add hline at the top
            .x_ref("x")
            .y_ref("paper")
            .line(ShapeLine::new()
                  .color(NamedColor::Black)
                  .width(0.7))
            );
    }

    fn gen_pband(selection: &Selection, cropped_projections: &Array5<f64>) -> Cube<f64> {
        let nspin    = cropped_projections.shape()[0];
        let nkpoints = cropped_projections.shape()[1];
        let nbands   = cropped_projections.shape()[2];

        let nspin = match nspin {
            1 | 4 => 1usize,
            2 => 2usize,
            _ => unreachable!(),
        };
        let mut summed_proj = Cube::<f64>::zeros([nspin, nkpoints, nbands]);

        if 1 == nspin {  // for ispin=1 or lsorbit=T system
            let proj = (0 .. nkpoints).into_par_iter()
                .map(|ik| {
                    let mut weights = Vector::zeros(nbands);
                    for ib in 0 .. nbands {
                        for is in selection.ispins.iter().cloned() {
                            for ia in selection.iatoms.iter().cloned() {
                                for iorbit in selection.iorbits.iter().cloned() {
                                    weights[ib] += cropped_projections[[is, ik, ib, ia, iorbit]];
                                }
                            }
                        }
                    }
                    weights
                })
                .collect::<Vec<Vector<f64>>>();

            for (ik, proj_item) in proj.iter().enumerate() {
                summed_proj.slice_mut(s![0usize, ik, ..]).assign(proj_item);
            }
        } else {  // ispin == 2
            for is in selection.ispins.iter().cloned() {
                let proj = (0 .. nkpoints).into_par_iter()
                    .map(|ik| {
                        let mut weights = Vector::zeros(nbands);
                        for ib in 0 .. nbands {
                            for ia in selection.iatoms.iter().cloned() {
                                for iorbit in selection.iorbits.iter().cloned() {
                                    weights[ib] += cropped_projections[[is, ik, ib, ia, iorbit]];
                                }
                            }
                        }
                        weights
                    })
                    .collect::<Vec<Vector<f64>>>();

                for (ik, proj_item) in proj.iter().enumerate() {
                    summed_proj.slice_mut(s![is, ik, ..]).assign(proj_item);
                }
            }
        };

        summed_proj
    }


    fn plot_pband(plot: &mut Plot, selection: &Selection, kpath: &Vector<f64>, cropped_eigvals: &Cube<f64>, projections: &Cube<f64>) {
        let nspin       = cropped_eigvals.shape()[0];
        let nkpoints    = cropped_eigvals.shape()[1];
        let nbands      = cropped_eigvals.shape()[2];

        assert_eq!(kpath.len(), nkpoints);      // cropped_eigvals[ispin, ikpoint, iband]

        let rand_color = RawSelection::get_random_color();
        let color = selection.color.clone().unwrap_or(rand_color.into());
        let marker = plotly::common::Marker::new().color(color);

        for ispin in 0 .. nspin {
            (0 .. nbands)
                .for_each(|iband| {
                    let dispersion = cropped_eigvals.slice(s![ispin, .., iband]).to_owned();
                    let projection = projections.slice(s![ispin, .., iband])
                        .iter()
                        .map(|x| {
                            if *x < 0.0 {
                                warn!("Negative projection number found: {} , it would be treated as zero", x);
                            }
                            (x * 20.0).ceil() as usize
                        })  // negative numbers are treated as 0
                        .collect::<Vec<usize>>();
                    let show_legend = 0 == iband && 0 == ispin;
                    let hover_template0 = match (nspin, ispin) {
                        (1, _) => format!("Band#: {}<br>", iband + 1),
                        (2, 0) => format!("Spin up<br>Band#: {}<br>", iband + 1),
                        (2, 1) => format!("Spin Down<br>Band#: {}<br>", iband + 1),
                        _ => unreachable!("Only two spin components available"),
                    };
                    let hover_template_array = projections.slice(s![ispin, .., iband])
                        .iter()
                        .map(|x| {
                            format!("{}E-Ef: %{{y:.4f}} eV<br>Projection: {:.3}", hover_template0, x)
                        })
                        .collect::<Vec<String>>();

                    let tr = Scatter::from_array(kpath.clone(), dispersion)
                        .mode(plotly::common::Mode::Markers)
                        .marker(marker.clone().opacity(0.4).size_array(projection))
                        .legend_group(&selection.label)
                        .show_legend(show_legend)
                        .hover_info(plotly::common::HoverInfo::Text)
                        .hover_template_array(hover_template_array)
                        .name(&selection.label);
                    plot.add_trace(tr);
                });
        }
    }


    fn gen_nclband(cropped_projections: &Array5<f64>, axis: Axis) -> Matrix<f64> {
        let nspin       = cropped_projections.shape()[0];
        let nkpoints    = cropped_projections.shape()[1];
        let nbands      = cropped_projections.shape()[2];

        assert_eq!(nspin, 4, "Not a NCL PROCAR");

        let iaxis = match axis {
            Axis::X => 1usize,
            Axis::Y => 2,
            Axis::Z => 3,
        };

        let mut summed_proj = Matrix::<f64>::zeros([nkpoints, nbands]);

        let proj = (0 .. nkpoints).into_par_iter()
            .map(|ik| {
                let mut weights = Vector::<f64>::zeros(nbands);
                for ib in 0 .. nbands {
                    weights[ib] += cropped_projections.slice(s![iaxis, ik, ib, .., ..]).sum();
                }
                weights
            })
            .collect::<Vec<Vector<f64>>>();

        for (ik, proj_item) in proj.iter().enumerate() {
            summed_proj.slice_mut(s![ik, ..]).assign(proj_item);
        }

        summed_proj
    }


    fn plot_nclband(plot: &mut Plot, kpath: &Vector<f64>, cropped_eigvals: &Cube<f64>, 
                     projections: &Matrix<f64>, colormap: plotly::common::ColorScalePalette, 
                     label: &str) {
        let nspin       = cropped_eigvals.shape()[0];
        let nkpoints    = cropped_eigvals.shape()[1];
        let nbands      = cropped_eigvals.shape()[2];

        assert_eq!(nspin, 1);
        assert_eq!(kpath.len(), nkpoints);

        (0 .. nbands)
            .for_each(|iband| {
                let dispersion = cropped_eigvals.slice(s![0, .., iband]).to_owned();
                let projection = projections.slice(s![.., iband]).to_owned().into_raw_vec_and_offset().0;
                let show_legend = 0 == iband;
                let hover_template_array = projection.iter()
                    .map(|x| {
                        format!("Band#: {}<Br>E-Ef: %{{y:.4f}} eV<br>{} Projection: {:.3}", iband + 1, label, x)
                    })
                    .collect::<Vec<String>>();

                let marker = plotly::common::Marker::new();
                /*
                 *let marker = if 0 == iband {
                 *    plotly::common::Marker::new()
                 *} else {
                 *    plotly::common::Marker::new()
                 *        //.color_bar(plotly::common::ColorBar::new()        // TODO: commented due to plotly-rs's stack overflow bug
                 *                   //.thickness(5)
                 *                   //.tick_vals(vec![-1.0, 1.0])
                 *                   //.outline_width(0))
                 *};
                 */

                let tr = Scatter::from_array(kpath.clone(), dispersion)
                    .mode(plotly::common::Mode::Markers)
                    .marker(marker
                            .color_scale(plotly::common::ColorScale::Palette(colormap.clone()))
                            .color_array(projection)
                            .cmin(-1.0)
                            .cmax(1.0))
                    .legend_group(label)
                    .show_legend(show_legend)
                    .hover_info(plotly::common::HoverInfo::Text)
                    .hover_template_array(hover_template_array)
                    .name(label);
                plot.add_trace(tr);
            });
    }


    // May be not useful here ...
    fn _filter_hse(procar: &mut Procar) -> bool {
        let skip_index = procar.kpoints.weights.iter()
            .position(|x| x.abs() < THRESHOLD);

        let skip_index = if let Some(i) = skip_index {
            i
        } else {
            return false;
        };


        procar.kpoints.nkpoints -= skip_index as u32;
        procar.kpoints.weights = procar.kpoints.weights
            .slice(s![skip_index ..])  // take weights[skip_index ..]
            .to_owned();
        procar.kpoints.kpointlist = procar.kpoints.kpointlist
            .slice(s![skip_index .., ..])
            .to_owned();

        procar.pdos.nkpoints = procar.kpoints.nkpoints;
        procar.pdos.eigvals = procar.pdos.eigvals
            .slice(s![.., skip_index .., ..])
            .to_owned();
        procar.pdos.occupations = procar.pdos.eigvals
            .slice(s![.., skip_index .., ..])
            .to_owned();
        procar.pdos.projected = procar.pdos.projected
            .slice(s![.., skip_index .., .., .., ..])
            .to_owned();

        let nkpoints = procar.kpoints.nkpoints as usize;

        assert!(
            procar.kpoints.weights.len()            == nkpoints &&
            procar.kpoints.kpointlist.shape()[0]    == nkpoints &&
            procar.pdos.nkpoints as usize           == nkpoints &&
            procar.pdos.eigvals.shape()[1]          == nkpoints &&
            procar.pdos.occupations.shape()[1]      == nkpoints &&
            procar.pdos.projected.shape()[1]        == nkpoints,
            "[*BUG*] Inconsistent k-point numbers in Procar instance"  // Treat as bug
            );

        true
    }
}


const TEMPLATE: &str = include_str!("./pband_template.toml");


impl OptProcess for Band {
    fn process(&self) -> Result<()> {
        if self.gen_template {
            let conf_filename = PathBuf::from("./pband.toml");

            info!("Generating projected band configuration tempalte ...");
            fs::write("pband.toml", TEMPLATE)?;
            info!("Template file written to {:?}. Exiting", &conf_filename);

            return Ok(());
        }

        
        // reading configuration
        let config = if let Some(config) = self.config.as_ref() {
            info!("Reading porjected band configuration fomr {:?}", self.config.as_ref());
            let config = fs::read_to_string(config)?;
            let config: Configuration = toml::from_str(&config)?;

            debug!("{:#?}", &config);

            Some(config)
        } else {
            None
        };

        let procar_fname    = config.as_ref().map(|cfg| &cfg.procar).unwrap_or(&self.procar);
        let outcar_fname    = config.as_ref().map(|cfg| &cfg.outcar).unwrap_or(&self.outcar);
        let txtout_prefix   = config.as_ref().map(|cfg| &cfg.txtout_prefix).unwrap_or(&self.txtout_prefix);
        let htmlout         = config.as_ref().map(|cfg| &cfg.htmlout).unwrap_or(&self.htmlout);
        let efermi          = config.as_ref().map(|cfg| &cfg.efermi).unwrap_or(&self.efermi);
        let ncl_spinor      = config.as_ref().map(|cfg| &cfg.ncl_spinor).unwrap_or(&self.ncl_spinor);
        let colormap        = config.as_ref().map(|cfg| &cfg.colormap).unwrap_or(&self.colormap);
        let kpoint_labels   = config.as_ref().map(|cfg| &cfg.kpoint_labels).unwrap_or(&self.kpoint_labels);
        let segment_ranges  = config.as_ref().map(|cfg| &cfg.segment_ranges).unwrap_or(&None);
        let ylim            = config.as_ref().map(|cfg| vec![cfg.ylim.0, cfg.ylim.1]).unwrap_or_else(|| self.ylim.clone());


        let mut procar: Result<Procar> = Err(anyhow!(""));
        let mut outcar: Result<Outcar> = Err(anyhow!(""));

        rayon::scope(|s| {
            s.spawn(|_| {
                info!("Reading band data from {:?}", procar_fname);
                procar = Procar::from_file(procar_fname);
            });
            s.spawn(|_| {
                info!("Reading fermi level and lattice data from {:?}", outcar_fname);
                outcar = Outcar::from_file(outcar_fname);
            });
        });

        let mut procar = procar.context(format!("Parse PROCAR file {:?} failed.", procar_fname))?;
        let outcar = outcar.context(format!("Parse OUTCAR file {:?} failed.", outcar_fname))?;


        let efermi = efermi.unwrap_or(outcar.efermi);
        let cell = outcar.ion_iters.last()
            .context("This OUTCAR doesn't complete at least one ionic step.")?
            .cell;

        let bcell = Poscar::acell_to_bcell(&cell).unwrap();

        info!("Found Fermi level: {}, shifting eigenvalues ...", efermi);
        procar.pdos.eigvals -= efermi;
        let procar = procar;  // rebind it, to remove mutability

        let nspin  = procar.pdos.nspin as usize;
        let is_ncl = procar.pdos.lsorbit;
        let nlm    = procar.pdos.nlm;
        let nions  = procar.pdos.nions as usize;
        let nbands = procar.pdos.nbands as usize;

        let segment_ranges      = segment_ranges.clone().unwrap_or(Self::find_segments(&procar.kpoints.kpointlist)?);
        let (kpath, kxs)        = Self::gen_kpath(&procar.kpoints.kpointlist, &bcell, &segment_ranges);
        let cropped_eigvals     = Self::gen_rawband(&procar.pdos.eigvals, &segment_ranges);
        let cropped_projections = Self::gen_cropped_projections(&procar.pdos.projected, &segment_ranges);

        let klabels = if let Some(label) = kpoint_labels.as_ref() {
            if label.len() != kxs.len() {
                bail!("Inconsistent k-point label number with segment ranges");
            }
            label.to_owned()
        } else {
            warn!("No k-point labels found, use empty labels instead");
            vec!["".to_string(); kxs.len()]
        };


        // Set up plot environment
        let mut plot = Plot::new();
        plot.use_local_plotly();

        let mut layout = plotly::Layout::new()
            .title(plotly::common::Title::with_text("Bandstructure"))
            .y_axis(plotly::layout::Axis::new()
                    .title(plotly::common::Title::with_text("E-Ef (eV)"))
                    .zero_line(true)
                    .range(ylim)
                    )
            .x_axis(plotly::layout::Axis::new()
                    .title(plotly::common::Title::with_text("Wavevector"))
                    .tick_values(kxs.clone())
                    .tick_text(klabels)
                    .zero_line(true)
                    )
            .height(960)
            .legend(plotly::layout::Legend::new().item_sizing(ItemSizing::Constant));

        Self::plot_boundaries(&mut layout, &kxs);
        plot.set_layout(layout);

        // Plot raw band
        info!("Plotting raw bands ...");
        Self::plot_rawband(&mut plot, kpath.clone(), &cropped_eigvals);

        // Plot ncl band
        if let Some(ax) = ncl_spinor.as_ref() {
            info!("Plotting ncl-band in {} direction", ax);
            let projected_band_ncl = Self::gen_nclband(&cropped_projections, *ax);
            let label = format!("Spinor {}", ax);
            Self::plot_nclband(&mut plot, &kpath, &cropped_eigvals, &projected_band_ncl, colormap.clone(), &label);

            let fname = PathBuf::from(&format!("{}_ncl_{}.txt", txtout_prefix, ax));
            let data = (0 .. nbands)
                .map(|iband| projected_band_ncl.slice(s![.., iband]).to_owned())
                .collect::<Vec<_>>();
            let data_ref = data.iter().collect::<Vec<&Vector<f64>>>();

            info!("Writing ncl band raw data to {:?} ...", &fname);
            write_array_to_txt(&fname, data_ref, "projection_coefficients nkpoints_x_nbands")?;
        }

        let selections = if config.as_ref().is_some() {
            if let Some(pband) = config.clone().unwrap().pband {
                Some(rawsel_to_sel(pband, nspin, is_ncl, &nlm, nions)?)
            } else {
                None
            }
        } else {
            None
        };

        // Plot projected bands
        if let Some(selections) = selections {
            info!("Plotting projected bands ...");
            let now = Instant::now();

            let pbands = selections
                .into_par_iter()
                .map(|sel| {
                    (
                        sel.clone(),
                        Self::gen_pband(&sel, &cropped_projections),
                    )
                })
                .collect::<Vec<_>>();
            
            for (sel, band) in pbands.into_iter() {
                Self::plot_pband(&mut plot, &sel, &kpath, &cropped_eigvals, &band);

                for is in &sel.ispins {
                    let spin_label = match (is_ncl, nspin, is) {
                        (false, 1, _) => {     "" },
                        (false, 2, 0) => {  "_up" },
                        (false, 2, 1) => {  "_dn" },
                        ( true, 1, 0) => { "_tot" },
                        ( true, 1, 1) => {  "_mx" },
                        ( true, 1, 2) => {  "_my" },
                        ( true, 1, 3) => {  "_mz" },
                        _ => { unreachable!("Invalied spin") },
                    };

                    let fname = PathBuf::from(&format!("{}_{}{}.txt", txtout_prefix, &sel.label, spin_label));
                    let data = (0 .. nbands)
                        .map(|iband| band.slice(s![*is, .., iband]).to_owned())
                        .collect::<Vec<_>>();
                    let data_ref = data.iter().collect::<Vec<&Vector<f64>>>();

                    info!("Writting projected band {} to {:?} ...", &sel.label, &fname);
                    write_array_to_txt(&fname, data_ref, "projection_coefficients nkpoints_x_nbands")?;
                }
            }

            info!("Projected band plot time usage: {:?}", now.elapsed());
        };


        // save data
        info!("Writing Bandstructure to {:?}", &htmlout);

        for is in 0 .. nspin {
            let spin_label = match (is_ncl, nspin, is) {
                (    _, 1, _) => { "" },
                (false, 2, 0) => { "_up" },
                (false, 2, 1) => { "_dn" },
                _ => { unreachable!("Invalied spin") },
            };

            let fname = PathBuf::from(&format!("{}{}.txt", txtout_prefix, spin_label));
            let mut data = (0 .. nbands)
                .map(|iband| cropped_eigvals.slice(s![is, .., iband]).to_owned())
                .collect::<Vec<_>>();
            data.insert(0, kpath.to_owned());

            let data_ref = data.iter().collect::<Vec<&Vector<f64>>>();
            info!("Writing raw band data to {:?}", &fname);
            write_array_to_txt(&fname, data_ref, "kpath(in_2pi) band-levels(nkpoints_x_nbands)")?;
        }

        plot.set_configuration(generate_plotly_configuration());
        plot.write_html(htmlout);

        if self.to_inline_html {
            info!("Printing inline html to stdout ...");
            println!("{}", plot.to_inline_html(None));
        }

        if self.show {
            plot.show();
        }
        
        Ok(())
    }
}


#[cfg(test)]
mod test {
    use super::*;
    use ndarray::{
        arr1,
        arr2,
    };

    const TEMPLATE_TEST: &str = include_str!("./pband_template_test.toml");

    #[test]
    fn test_gen_kpath() {
        let kpoints = arr2(&[[ 0.     ,  0.     ,  0.     ],  // 1
                             [-0.08325,  0.16675,  0.     ],  // 2
                             [-0.1665 ,  0.3335 ,  0.     ],  // 3
                             [-0.24975,  0.50025,  0.     ],  // 4
                             [-0.333  ,  0.667  ,  0.     ],  // 5 <- sep
                             [-0.333  ,  0.667  ,  0.     ],  // 6 <- sep
                             [-0.24975,  0.62525,  0.     ],  // 7
                             [-0.1665 ,  0.5835 ,  0.     ],  // 8
                             [-0.08325,  0.54175,  0.     ],  // 9
                             [ 0.     ,  0.5    ,  0.     ],  // 10 <- sep
                             [ 0.     ,  0.5    ,  0.     ],  // 11 <- sep
                             [ 0.     ,  0.375  ,  0.     ],  // 12
                             [ 0.     ,  0.25   ,  0.     ],  // 13
                             [ 0.     ,  0.125  ,  0.     ],  // 14
                             [ 0.     ,  0.     ,  0.     ]]);// 15

        let cell = [[    2.884_730_620_005_97 * 0.993,  -1.665_5 * 0.993,   0.0000000000000000 * 0.993],
                    [    0.0000000000000000 * 0.993,   3.331 * 0.993,   0.0000000000000000 * 0.993],
                    [    0.0000000000000000 * 0.993,   0.0000000000000000 * 0.993,  23.164 * 0.993]];

        let bcell = Poscar::acell_to_bcell(&cell).unwrap();

        let (kpath, _) = Band::gen_kpath(&kpoints, &bcell, &[(1, 5), (6, 10), (11, 15), (6, 10)]);
        let expect = arr1(&[0.0, 0.05041295144327735, 0.1008259028865547, 0.15123885432983203, 0.2016518057731094,
                          0.2016518057731094, 0.22682051907635853, 0.2519892323796077, 0.27715794568285684, 0.30232665898610595,
                          0.30232665898610595, 0.3459637207291467, 0.38960078247218755, 0.4332378442152283, 0.4768749059582691,
                          0.4768749059582691, 0.5020436192615182, 0.5272123325647673, 0.5523810458680164, 0.5775497591712655]);
        eprintln!("{}", &kpath);
        assert!((kpath - &expect).iter().all(|x| x.abs() < 1E-6));

        let (kpath, _) = Band::gen_kpath(&kpoints, &bcell, &[(5, 1), (10, 6), (15, 11), (10, 6)]);
        assert!((kpath - &expect).iter().all(|x| x.abs() < 1E-6));
        

        let segment_ranges = Band::find_segments(&kpoints).unwrap();
        assert_eq!(segment_ranges, vec![(1, 5), (6, 10), (11, 15)]);
    }


    #[test]
    fn test_config() {
        let nlm = "s     py     pz     px    dxy    dyz    dz2    dxz    dx2" 
            .split_whitespace()
            .map(String::from)
            .collect::<Vec<_>>();
        let nspin = 2;
        let nions = 8usize;
        let is_ncl = false;
        let kpoint_labels_ref = vec!["G", "K", "M", "G"].into_iter()
            .map(String::from)
            .collect::<Vec<String>>();

        let c: Configuration = toml::from_str(TEMPLATE_TEST).unwrap();
        let v = rawsel_to_sel(c.clone().pband.unwrap(), nspin, is_ncl, &nlm, nions).unwrap();

        assert_eq!(c.kpoint_labels.as_ref(), Some(&kpoint_labels_ref));
        assert_eq!(c.txtout_prefix, "band_raw");
        assert_eq!(c.segment_ranges, Some(vec![(1, 40), (41, 80), (81, 120)]));
        assert_eq!(c.ncl_spinor.unwrap(), Axis::Z);
        //assert_eq!(c.colormap, ColorScalePalette::Blues);  // Palette doesn't derive Eq trait
        assert_eq!(c.efermi.as_ref(), Some(&0.0));

        assert_eq!(v[0].label, "plot1");
        assert_eq!(v[0].iatoms, &[0, 2, 3, 4, 5, 6, 7]);
        assert_eq!(v[0].iorbits, &[0, 3, 4]);
        assert_eq!(c.ylim, (-2.0, 5.0));

        let s = toml::to_string(&c).unwrap();
        println!("{}", s);
        println!("{:?}", v);
    }
}
