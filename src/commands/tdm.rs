use std::{
    fs,
    fmt::Write,
    path::PathBuf,
};

use clap::{
    Parser,
    AppSettings,
};
use anyhow::bail;
use log::{
    warn,
    info,
};
use itertools::iproduct;
use rayon::prelude::*;
use plotly;

use crate::{
    types::{
        Result,
        Axis,
        Vector,
    },
    OptProcess,
    vasp_parsers::wavecar::{
        Wavecar,
        WavecarType,
    },
    commands::common::write_array_to_txt,
};


#[derive(Debug, Parser)]
#[clap(setting = AppSettings::ColoredHelp,
       setting = AppSettings::ColorAuto)]
/// Calculate Transition Dipole Moment (TDM) between given bands.
///
/// Note: This command can only calculate the TDM between bands in
/// same k-point. The inter-kpoint transition is not supported yet.
/// Also, this commands calculates the TDM in reciprocal space by
///
/// tdm_{i->j} = <phi_j|e*r|phi_i> = i*ħ/(ΔE*m)*<phi_j|p|phi_i>
pub struct Tdm {
    #[clap(short, long, default_value = "./WAVECAR")]
    /// WAVECAR file path.
    wavecar: PathBuf,

    #[clap(long, possible_values = &["x", "z"])]
    /// Gamma Half direction of WAVECAR. You need to set this to 'x' or 'z' when
    /// processing WAVECAR produced by `vasp_gam`.
    gamma_half: Option<String>,

    #[clap(short = 's', long, default_value = "1", possible_values = &["1", "2"])]
    /// Spin index, 1 for up, 2 for down.
    ispin: usize,

    #[clap(short = 'k', long, default_value = "1")]
    /// K-point index, starts from 1.
    ikpoint: usize,

    #[clap(short = 'i', long, min_values = 1, required = true)]
    /// Initial band indices, start from 1.
    ibands: Vec<usize>,

    #[clap(short = 'j', long, min_values = 1, required = true)]
    /// Final band indices, starts from 1.
    jbands: Vec<usize>,

    #[clap(long, default_value = "0.05")]
    /// Smearing width, in eV.
    sigma: f64,

    #[clap(short, long)]
    /// Print the calculated TDM to screen.
    verbose: bool,

    #[clap(long, default_value = "tdm_peaks.txt")]
    /// Write the TDM peaks to raw txt file.
    peakout: PathBuf,

    #[clap(long, default_value = "tdm_smeared.txt")]
    /// Write the summed and smeared TDM to raw txt file.
    txtout: PathBuf,

    #[clap(long, default_value = "tdm_smeared.html")]
    /// Write the plot of TDM to html file.
    htmlout: PathBuf,

    #[clap(long)]
    /// Print the inline HTML to stdout.
    to_inline_html: bool,

    #[clap(long)]
    /// Open the default browser to show the plot.
    show: bool,

    #[clap(long, default_value = "0.1")]
    /// Specify the width of bars in the center of peaks. (eV)
    barwidth: f64,
}


impl Tdm {
    fn check_and_transform_band_index(ibands: &[usize], nbands: usize) -> Result<Vec<usize>> {
        let mut ibands = ibands.to_vec();
        ibands.sort_unstable();
        ibands.dedup();

        for iband in ibands.iter().cloned() {
            if iband == 0 || iband > nbands {
                bail!("Invalid band index: {} . It should be >= 1 and <= {}.", iband, nbands);
            }
        }

        // indices start from 0
        Ok(ibands.into_iter().map(|i| i - 1).collect())
    }


    // lorentz_smearing(x::AbstractArray, x0::Float64, Γ=0.05) = @. Γ/(2π) / ((x-x0)^2 + (Γ/2)^2)
    fn smearing_lorentz(x: &[f64], x0s: &[f64], gamma: f64, scales: &[f64]) -> Vector<f64> {
        const PI: f64 = std::f64::consts::PI;

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

    fn apply_smearing(x: &[f64], centers: &[f64], width: f64, scales: Option<&[f64]>) -> Vector<f64> {
        let clen = centers.len();
        let mut fac = vec![1.0; 0];

        let scales = if let Some(factors) = scales {
            assert_eq!(factors.len(), clen, "[DOSsmear]: factors' length inconsistent with smearing peaks");
            factors
        } else {
            fac.resize(clen, 1.0);
            &fac
        };

        Self::smearing_lorentz(x, centers, width, scales)
    }

    fn gen_smeared_tdm(x: &[f64], tdms: Vec<(usize, usize, f64, f64, f64, f64, f64, f64)>, sigma: f64) -> Vec<Vector<f64>> {
        let mut smeared_tdms = vec![];  // x, y, z, tot

        let centers = tdms.iter().map(|t| t.4).collect::<Vec<f64>>();
        let txs     = tdms.iter().map(|t| t.5).collect::<Vec<f64>>();
        let tys     = tdms.iter().map(|t| t.6).collect::<Vec<f64>>();
        let tzs     = tdms.iter().map(|t| t.7).collect::<Vec<f64>>();

        smeared_tdms.push(Self::apply_smearing(&x, &centers, sigma, Some(&txs)));
        smeared_tdms.push(Self::apply_smearing(&x, &centers, sigma, Some(&tys)));
        smeared_tdms.push(Self::apply_smearing(&x, &centers, sigma, Some(&tzs)));
        let tot_tdms = smeared_tdms[0].clone() + &smeared_tdms[1] + &smeared_tdms[2];
        smeared_tdms.push(tot_tdms);

        smeared_tdms
    }

}




impl OptProcess for Tdm {
    fn process(&self) -> Result<()> {
        info!("Reading {:?}", &self.wavecar);
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


        let ispin = if self.ispin > 0 && self.ispin <= wav.nspin as usize {
            self.ispin - 1
        } else {
            bail!("Invalid ispin: ispin shoule be >= 1 and <= {}.", wav.nspin);
        };

        let ikpoint = if self.ikpoint > 0 && self.ikpoint <= wav.nkpoints as usize {
            self.ikpoint - 1
        } else {
            bail!("Invalid ikpoint: ikpoint should be >= 1 and <= {}.", wav.nkpoints);
        };


        let nbands = wav.nbands as usize;
        let eigs   = &wav.band_eigs;
        let efermi = wav.efermi;
        let ibands = Self::check_and_transform_band_index(self.ibands.as_slice(), nbands)?;
        let jbands = Self::check_and_transform_band_index(self.jbands.as_slice(), nbands)?;

        let tdms = iproduct!(ibands, jbands)
            .filter(|(iband, jband)| iband < jband)
            .par_bridge()
            .map(|(iband, jband)| {
                let eig_i = eigs[(ispin, ikpoint, iband)] - efermi;
                let eig_j = eigs[(ispin, ikpoint, jband)] - efermi;
                #[allow(non_snake_case)]
                let dE    = eig_j - eig_i;

                let [tdmx, tdmy, tdmz] = wav.transition_dipole_moment(ispin as u64, ikpoint as u64, iband as u64, jband as u64);

                let tdmx = tdmx.norm();
                let tdmy = tdmy.norm();
                let tdmz = tdmz.norm();

                (iband, jband, eig_i, eig_j, dE, tdmx, tdmy, tdmz)
            })
            .collect::<Vec<_>>();

        {
            let mut txt: String = String::new();
            
            writeln!(&mut txt, "# iband jband     E_i     E_j      ΔE        Tx       Ty       Tz")?;
            tdms.iter()
                .map(|(iband, jband, eig_i, eig_j, d_e, tx, ty, tz)| {
                    writeln!(&mut txt, "  {:5} {:5} {:7.3} {:7.3} {:7.3}  {:8.3} {:8.3} {:8.3}",
                             iband+1, jband+1, eig_i, eig_j, d_e, tx, ty, tz)?;
                    Ok(())
                })
                .collect::<Result<()>>()?;

            if self.verbose {
                print!("{}", txt);
            }

            info!("Writing peaks info to {:?} ...", &self.peakout);
            fs::write(&self.peakout, &txt)?;
        }


        // Plot with plotly
        let mut plot = plotly::Plot::new();

        let des = tdms.iter().map(|t| t.4).collect::<Vec<f64>>();
        let txs = tdms.iter().map(|t| t.5).collect::<Vec<f64>>();
        let tys = tdms.iter().map(|t| t.6).collect::<Vec<f64>>();
        let tzs = tdms.iter().map(|t| t.7).collect::<Vec<f64>>();

        let hover_template_array = tdms.iter()
            .map(|(iband, jband, _, _, de, _, _, _)| {
                format!("{}->{},dE={:.3}eV<br>%{{y:.4f}}Debye", iband+1, jband+1, de)
            })
            .collect::<Vec<String>>();

        for (label, t, color) in [("Tx", txs, "#1f77b4"), ("Ty", tys, "#ff7f0e"), ("Tz", tzs, "#2ca02c")] {
            let tr = plotly::Bar::new(des.clone(), t)
                .marker(plotly::common::Marker::new().color(color))
                .width(self.barwidth)
                .name(label)
                .hover_template_array(hover_template_array.clone())
                .legend_group(label);
            plot.add_trace(tr);
        }


        let x_min = tdms.iter().map(|t| t.4).reduce(f64::min).unwrap() - 2.0;
        let x_max = tdms.iter().map(|t| t.4).reduce(f64::max).unwrap() + 2.0;
        let nx = (x_max - x_min).ceil() as usize * 500;
        let x = Vector::<f64>::linspace(x_min, x_max, nx);
        let smeared_tdms = Self::gen_smeared_tdm(&x.to_vec(), tdms, self.sigma);
        
        // Write smeared TDM to txt
        {
            let dat = std::iter::once(&x).chain(smeared_tdms.iter()).collect::<Vec<_>>();
            info!("Writing smeared TDM data to {:?}", &self.txtout);
            write_array_to_txt(&self.txtout, dat, "E(eV)    Tx(Debye)   Ty(Debye)   Tz(Debye)")?;
        }

        for (label, i, color) in [("Tx", 0, "#1f77b4"), ("Ty", 1, "#ff7f0e"), ("Tz", 2, "#2ca02c"), ("T", 3, "000000")] {
            let tr = plotly::Scatter::from_array(x.clone(), smeared_tdms[i].clone())
                .mode(plotly::common::Mode::Lines)
                .hover_info(plotly::common::HoverInfo::None)
                .name(label)
                .legend_group(label)
                .marker(plotly::common::Marker::new().color(color));

            plot.add_trace(tr);
        }

        plot.use_local_plotly();
        let layout = plotly::Layout::new()
            .bar_mode(plotly::layout::BarMode::Stack)
            .title(plotly::common::Title::new("Transition Dipole Moments"))
            .y_axis(plotly::layout::Axis::new()
                    .title(plotly::common::Title::new("TDM (Debye)"))
                    .fixed_range(false))
            .x_axis(plotly::layout::Axis::new()
                    .title(plotly::common::Title::new("Energy (eV)"))
                    .fixed_range(false))
            .hover_mode(plotly::layout::HoverMode::X);

        plot.set_layout(layout);

        // Write html
        info!("Writing plot to {:?}", &self.htmlout);
        plot.to_html(&self.htmlout);

        if self.to_inline_html {
            info!("Printing inline HTML to stdout ...");
            println!("{}", plot.to_inline_html(None));
        }

        if self.show {
            plot.show();
        }

        Ok(())
    }
}
