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
/// Plot wavefunction in realspace, and save it as '.vasp' file.
pub struct Wav3D {
    #[structopt(long, short="w", default_value="./WAVECAR")]
    /// WAVECAR file name.
    wavecar: PathBuf,

    #[structopt(long, short="p", default_value="./POSCAR")]
    /// POSCAR filename, POSCAR is needed to get the real-space wavefunction.
    poscar: PathBuf,

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

    #[structopt(long, short="o", possible_values=&["normsquared", "ns", "real", "re", "imag", "im", "reim"])]
    /// Specify output part of the wavefunction.
    ///
    /// Detailed message:{n}
    /// - normsquared/ns: Perform `ρ(r) = |ѱ(r)|^2` action to get the spatial distribution of selected band.{n}
    /// - real/re: Real part of the wavefunction, suffix '_re.vasp' is added to the output filename.{n}
    /// - imag/im: Imaginary part of the wavefunction, suffix '_im.vasp' is added to the output filename.{n}
    /// - reim: Output both real part and imaginary parts of the wavefunction.
    output_parts: Vec<String>,

    #[structopt(long, default_value="wav")]
    /// Prefix of output filename.
    prefix: String,

    #[structopt(long, short="e")]
    /// Add eigen value suffix to the filename
    show_eigs_suffix: bool,
}


fn save_to_vasp(fname: &str, chgd: ndarray::Array3<f64>, pos: &Poscar) -> Result<()> {
    let ngrid = [chgd.raw_dim()[0], chgd.raw_dim()[1], chgd.raw_dim()[2]];

    let fname = PathBuf::from(fname);
    if fname.is_file() {
        warn!("File {:?} exists, overwriting ...", fname);
    } else {
        info!("Writing {:?} ...", fname);
    }
    chg::ChargeDensity {
        chgtype: chg::ChargeType::Locpot,
        pos: pos.clone(),
        ngrid,
        chg: vec![chgd],
        aug: vec![],
    }.to_file(&fname)
}


impl OptProcess for Wav3D {
    fn process(&self) -> Result<()> {
        println!("{:?}", self);

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

        info!("Reading POSCAR: {:?}", &self.poscar);
        let pos = Poscar::from_file(&self.poscar)?;

        let ngrid = wav.ngrid;
        let efermi = wav.efermi;
        let eigs   = wav.band_eigs.clone();
        let factor = ngrid.iter().product::<u64>() as f64 * 8.0;

        let has_normsquared = self.output_parts.iter().any(|s| s == "normsquared" || s == "ns");
        let has_real = self.output_parts.iter().any(|s| s == "real" || s == "re" || s == "reim");
        let has_imag = self.output_parts.iter().any(|s| s == "imag" || s == "im" || s == "reim");

        if !(has_normsquared || has_real || has_imag) {
            warn!("You have not specify the `output_parts` or `list`, rsgrad did nothing.");
            return Ok(())
        }


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
        let wav = wav;

        indices.into_par_iter()
            .map(|(ispin, ikpoint, iband)| {
                info!("Processing spin {}, k-point {:3}, band {:4} ...", ispin, ikpoint, iband);

                let eigs_suffix = if self.show_eigs_suffix {
                    format!("_{:06.3}eV", eigs[[ispin as usize, ikpoint as usize, iband as usize]] - efermi)
                } else {
                    String::new()
                };

                let wavr = wav.get_wavefunction_realspace(ispin, ikpoint, iband)?.normalize();
                let chgd = match wavr.clone() {
                    Wavefunction::Complex64Array3(w)  => w.mapv(|v| v.norm_sqr() * factor),
                    Wavefunction::Float64Array3(w)    => w.mapv(|v| v * v * factor),
                    Wavefunction::Ncl64Array4(w)      => {
                        w.slice(s![0usize, .., .., ..]).mapv(|v| v.norm_sqr() * factor) +
                            w.slice(s![1usize, .., .., ..]).mapv(|v| v.norm_sqr() * factor)
                    },
                    _ => unreachable!("Invalid Wavefunction type."),
                };

                if has_normsquared {
                    let ofname = format!("{}_{}-{}-{}{}.vasp", &self.prefix, ispin+1, ikpoint+1, iband+1, eigs_suffix);
                    save_to_vasp(&ofname, chgd, &pos)?;
                }

                let (d1, d2, d3, d4) = match wavr {
                    Wavefunction::Complex64Array3(w) => ( Some(w.mapv(|x| x.re * factor)), Some(w.mapv(|x| x.im * factor)), None, None ),
                    Wavefunction::Float64Array3(w)   => ( Some(w), None, None, None),
                    Wavefunction::Ncl64Array4(w)     => (
                        Some(w.slice(s![0usize, .. ,.. ,..]).mapv(|v| v.re * factor)),
                        Some(w.slice(s![0usize, .. ,.. ,..]).mapv(|v| v.im * factor)),
                        Some(w.slice(s![1usize, .. ,.. ,..]).mapv(|v| v.re * factor)),
                        Some(w.slice(s![1usize, .. ,.. ,..]).mapv(|v| v.im * factor)),
                        ),
                    _ => unreachable!("Invalid Wavefunction type."),
                };

                if has_real {
                    match wavecar_type {
                        WavecarType::NonCollinear => {
                            let ofname = format!("{}_{}-{}-{}{}_ure.vasp", &self.prefix, ispin+1, ikpoint+1, iband+1, eigs_suffix);
                            save_to_vasp(&ofname, d1.unwrap(), &pos)?;
                            let ofname = format!("{}_{}-{}-{}{}_dre.vasp", &self.prefix, ispin+1, ikpoint+1, iband+1, eigs_suffix);
                            save_to_vasp(&ofname, d3.unwrap(), &pos)?;
                        },
                        _ => {
                            let ofname = format!("{}_{}-{}-{}{}_re.vasp", &self.prefix, ispin+1, ikpoint+1, iband+1, eigs_suffix);
                            save_to_vasp(&ofname, d1.unwrap(), &pos)?;
                        },
                    }
                }

                if has_imag {
                    match wavecar_type {
                        WavecarType::Standard => {
                            let ofname = format!("{}_{}-{}-{}{}_im.vasp", &self.prefix, ispin+1, ikpoint+1, iband+1, eigs_suffix);
                            save_to_vasp(&ofname, d2.unwrap(), &pos)?;
                        },
                        WavecarType::GamaHalf(_) => {
                            bail!("Gamma-halved wavefunction doesn't have imaginary part, please check your input.");
                        },
                        WavecarType::NonCollinear => {
                            let ofname = format!("{}_{}-{}-{}{}_uim.vasp", &self.prefix, ispin+1, ikpoint+1, iband+1, eigs_suffix);
                            save_to_vasp(&ofname, d2.unwrap(), &pos)?;
                            let ofname = format!("{}_{}-{}-{}{}_dim.vasp", &self.prefix, ispin+1, ikpoint+1, iband+1, eigs_suffix);
                            save_to_vasp(&ofname, d4.unwrap(), &pos)?;
                        },
                    }
                }

                Ok(())
            })
            .collect::<Result<()>>()
    }
}
