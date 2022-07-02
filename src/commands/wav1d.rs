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

use crate::{
    vasp_parsers::{
        chg,
        wavecar::{
            Wavecar,
            WavecarType,
            Wavefunction,
        },
    },
    types::{
        Result,
        OptProcess,
        Axis,
    },
    Poscar,
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
    /// Specify the file to be written with raw wav1d data.
    txtout: PathBuf,

    #[structopt(long, default_value="z",
                possible_values = &Axis::variants(),
                case_insensitive = true)]
    /// Integration direction. e.g. if 'z' is provided, the XoY plane is integrated.
    axis: Axis,

    #[structopt(long)]
    /// Render the plot and print thw rendered code to stdout.
    to_inline_html: bool,
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
        } else {
            if wav.wavecar_type != WavecarType::Standard &&
               wav.wavecar_type != WavecarType::NonCollinear {
                warn!("Current WAVECAR is gamma-halved, sometimes the gamma-x and gamma-z verions have same plane wave numbers.
I suggest you provide `gamma_half` argument to avoid confusion.");
            }
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

        let wavecar_type = wav.wavecar_type.clone();
        let wav = wav;  // Cancel the mutability

        let mut dat = indices.into_par_iter()
            .map(|(ispin, ikpoint, iband)| {
                info!("Processing spin {}, k-point {:3}, band {:4} ...", ispin+1, ikpoint+1, iband+1);
                let eig = eigs[[ispin as usize, ikpoint as usize, iband as usize]] - efermi;
                let label = format!("s{} k{} b{} {:06.3}eV", ispin+1, ikpoint+1, iband+1, eig);

                let wavr = wav.get_wavefunction_realspace(ispin, ikpoint, iband)
                    .expect(&format!("Failed to get wavefunction in realspace at s{} k{} b{}", ispin+1, ikpoint+1, iband+1))
                    .normalize();

                let chgd = match wavr.clone() {
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
                        chgd.mean_axis(ndarray::Axis(2)).unwrap()
                            .mean_axis(ndarray::Axis(1)).unwrap()
                    },
                    Axis::Y => {
                        chgd.mean_axis(ndarray::Axis(2)).unwrap()
                            .mean_axis(ndarray::Axis(0)).unwrap()
                    },
                    Axis::Z => {
                        chgd.mean_axis(ndarray::Axis(1)).unwrap()
                            .mean_axis(ndarray::Axis(0)).unwrap()
                    },
                };

                (eig, label, chg1d)
            })
            .collect::<Vec<(f64, String, Array1<f64>)>>();

        //dat.sort_unstable_by(|()|)



        Ok(())
    }
}
