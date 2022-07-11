use std::path::PathBuf;

use structopt::{
    StructOpt,
    clap::AppSettings,
};
use log::{
    info,
    warn,
};
use anyhow::bail;
use ndarray::s;
use rayon::prelude::*;
use itertools::iproduct;
use ndarray::Array1;
use plotly;

use crate::{
    vasp_parsers::wavecar::{
        Wavecar,
        WavecarType,
        Wavefunction,
    },
    types::{
        Result,
        OptProcess,
        Axis,
    },
    commands::common::write_array_to_txt,
};


#[derive(Debug, StructOpt)]
#[structopt(setting = AppSettings::ColoredHelp,
            setting = AppSettings::ColorAuto,
            setting = AppSettings::AllowNegativeNumbers)]
/// Plot wavefunction in realspace, then integrate over some plane, and save it as '.txt' file.
pub struct Wav1D {
    #[structopt(long, short="w", default_value="./WAVECAR")]
    /// WAVECAR file name.
    wavecar: PathBuf,

    #[structopt(long, short="s", default_value="1")]
    /// Select spin index, starting from 1.
    ispins: Vec<i32>,

    #[structopt(long, short="k", default_value="1")]
    /// Select kpoint index, starting from 1.
    ikpoints: Vec<i32>,

    #[structopt(long, short="b")]
    /// Select band index, starting from 1.
    ibands: Vec<i32>,

    #[structopt(long, short="l")]
    /// List the brief info of current WAVECAR.
    list: bool,

    #[structopt(long, possible_values=&["x", "z"])]
    /// Gamma Half direction of WAVECAR. You need to set this to 'x' or 'z' when
    /// processing WAVECAR produced by `vasp_gam`.
    gamma_half: Option<String>,

//  #[structopt(long, short="o", default_value="ns",
//              possible_values=&["normsquared", "ns", "real", "re", "imag", "im"])]
//  /// Specify output part of the wavefunction.
//  ///
//  /// Detailed message:{n}
//  /// - normsquared/ns: Perform `ρ(r) = |ѱ(r)|^2` action to get the spatial distribution of selected band.{n}
//  /// - real/re: Real part of the wavefunction, suffix '_re.vasp' is added to the output filename.{n}
//  /// - imag/im: Imaginary part of the wavefunction, suffix '_im.vasp' is added to the output filename.{n}
//  /// - reim: Output both real part and imaginary parts of the wavefunction.
//  output_part: String,

    #[structopt(long, default_value="wav1d.txt")]
    /// Specify the file name to be written with raw wav1d data.
    txtout: PathBuf,

    #[structopt(long, default_value="wav1d.html")]
    /// Specify the file name to be written with html wav1d data.
    htmlout: PathBuf,

    #[structopt(long, default_value="z",
                possible_values = &Axis::variants(),
                case_insensitive = true)]
    /// Integration direction. e.g. if 'z' is provided, the XoY plane is integrated.
    axis: Axis,

    #[structopt(long)]
    /// Render the plot and print thw rendered code to stdout.
    to_inline_html: bool,

    #[structopt(long)]
    /// Open the browser and show the plot immediately.
    show: bool,

    #[structopt(long, default_value="10")]
    /// Scale the wavefunction.
    scale: f64,
}


impl OptProcess for Wav1D {
    fn process(&self) -> Result<()> {
        info!("Reading WAVECAR: {:?}", &self.wavecar);
        let mut wav = Wavecar::from_file(&self.wavecar)?;
        if let Some(gammahalf) = self.gamma_half.as_ref() {
            if wav.wavecar_type == WavecarType::Standard ||
               wav.wavecar_type == WavecarType::NonCollinear {
                    bail!("Current WAVECAR is not gamma-halved, rsgrad can determine the WAVECAR type directly, \
please remove the argument `gamma_half`.")
            }

            let gammahalf = match gammahalf.as_ref() {
                "x" => WavecarType::GamaHalf(Axis::X),
                "z" => WavecarType::GamaHalf(Axis::Z),
                _ => panic!("Unreachable branch"),
            };
            
            wav.set_wavecar_type(gammahalf)?;
        } else if wav.wavecar_type != WavecarType::Standard &&
               wav.wavecar_type != WavecarType::NonCollinear {
                warn!("Current WAVECAR is gamma-halved, sometimes the gamma-x and gamma-z verions have same plane wave numbers.
I suggest you provide `gamma_half` argument to avoid confusion.");
        }

        if self.list {
            println!("{}", wav);
        }

        let efermi = wav.efermi;
        let eigs   = wav.band_eigs.clone();

        let ispins = self.ispins.iter()
            .cloned()
            .map(|v| v as u64 - 1)
            .collect::<Vec<_>>();
        let ikpoints = self.ikpoints.iter()
            .cloned()
            .map(|v| v as u64 - 1)
            .collect::<Vec<_>>();
        let ibands = self.ibands.iter()
            .cloned()
            .map(|v| v as u64 - 1)
            .collect::<Vec<_>>();

        let indices = iproduct!(ispins, ikpoints, ibands)
            .collect::<Vec<(u64, u64, u64)>>();

        let wav = wav;  // Cancel the mutability

        let mut dat = indices.into_par_iter()
            .map(|(ispin, ikpoint, iband)| {
                info!("Processing spin {}, k-point {:3}, band {:4} ...", ispin+1, ikpoint+1, iband+1);
                let eig = eigs[[ispin as usize, ikpoint as usize, iband as usize]] - efermi;
                let label = format!("s{}_k{}_b{}_{:06.3}eV", ispin+1, ikpoint+1, iband+1, eig);

                let wavr = wav.get_wavefunction_realspace(ispin, ikpoint, iband)
                    .unwrap_or_else(|_| panic!("Failed to get wavefunction in realspace at s{} k{} b{}", ispin+1, ikpoint+1, iband+1))
                    .normalize();

                let chgd = match wavr {
                    Wavefunction::Complex64Array3(w)  => w.mapv(|v| v.norm_sqr()),
                    Wavefunction::Float64Array3(w)    => w.mapv(|v| v * v),
                    Wavefunction::Ncl64Array4(w)      => {
                        w.slice(s![0usize, .., .., ..]).mapv(|v| v.norm_sqr()) +
                            w.slice(s![1usize, .., .., ..]).mapv(|v| v.norm_sqr())
                    },
                    _ => unreachable!("Invalid Wavefunction type."),
                };

                let chg1d = match self.axis {
                    Axis::X => {
                        chgd.sum_axis(ndarray::Axis(2))
                            .sum_axis(ndarray::Axis(1))
                    },
                    Axis::Y => {
                        chgd.sum_axis(ndarray::Axis(2))
                            .sum_axis(ndarray::Axis(0))
                    },
                    Axis::Z => {
                        chgd.sum_axis(ndarray::Axis(1))
                            .sum_axis(ndarray::Axis(0))
                    },
                } * self.scale;

                (eig, label, chg1d)
            })
            .collect::<Vec<(f64, String, Array1<f64>)>>();

        // sort in descending order
        dat.sort_unstable_by(|(e1, _, _), (e2, _, _)| e2.partial_cmp(e1).unwrap());
        
        let iaxis = match self.axis {
            Axis::X => 0usize,
            Axis::Y => 1usize,
            Axis::Z => 2usize,
        };
        let axislen = {
            let r = wav.acell[iaxis];
            (r[0] + r[0] + r[1] * r[1] + r[2] * r[2]).sqrt()
        };
        let xdat = ndarray::Array::linspace(0.0, axislen, wav.ngrid[iaxis] as usize * 2);
        
        let mut plot = plotly::Plot::new();

        dat.iter()
            .for_each(|(e, l, w)| {
                let trace = plotly::Scatter::from_array(xdat.clone(), w.mapv(|x| x+e))
                    .mode(plotly::common::Mode::Lines)
                    .name(l);
                plot.add_trace(trace);
            });

        let layout = plotly::Layout::new()
            .title(plotly::common::Title::new(&format!("Wavefunction Along {} Axis", self.axis)))
            .y_axis(plotly::layout::Axis::new()
                    .title(plotly::common::Title::new("E-E<sub>f</sub> (eV)"))
                    .zero_line(true))
            .x_axis(plotly::layout::Axis::new()
                    .title(plotly::common::Title::new("Distance (Å)")));
        plot.set_layout(layout);

        info!("Writing to {:?}", self.htmlout);
        plot.to_html(&self.htmlout);

        let comment = dat.iter()
            .map(|(_, l, _)| l.clone())
            .collect::<Vec<String>>()
            .join(" ");
        let header = format!("Distance(A) {}", comment);
        let data_ref = std::iter::once(&xdat)
            .chain(dat.iter().map(|(_, _, w)| w))
            .collect::<Vec<&Array1<f64>>>();

        info!("Writing to {:?}", self.txtout);
        write_array_to_txt(&self.txtout, data_ref, &header)?;

        if self.show {
            plot.show();
        }

        if self.to_inline_html {
            info!("Printing inline html to stdout ...");
            println!("{}", plot.to_inline_html(None));
        }


        Ok(())
    }
}
