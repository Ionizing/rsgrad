use std::path::PathBuf;

use clap::{
    Parser,
    AppSettings,
};
use log::{
    info,
    warn,
};
use ndarray::{
    Array1,
    Array2,
    Array3,
    Axis,
};
use anyhow::Context;
use itertools::multizip;
use colored::Colorize;

use crate::{
    vasp_parsers::{
        wavecar::Wavecar,
        procar::Procar,
    },
    types::{
        Result,
        OptProcess,
    },
};


#[derive(Debug, Parser)]
#[clap(setting = AppSettings::ColoredHelp,
       setting = AppSettings::ColorAuto)]
/// Find band gap and print positions of VBM and CBM
pub struct Gap {
    #[clap(long, short = 'w', default_value = "WAVECAR")]
    /// WAVECAR file name
    wavecar: PathBuf,

    #[clap(long, short = 'p', default_value = "PROCAR")]
    /// PROCAR file name
    procar: PathBuf,
}


impl Gap {
    fn bands_from_procar(procar: &Procar) -> (Array3<f64>, Array3<f64>, Array2<f64>) {
        let eigs = procar.pdos.eigvals.clone();
        let occs = procar.pdos.occupations.clone();
        let kvec = procar.kpoints.kpointlist.clone();

        (eigs, occs, kvec)
    }

    fn bands_from_wavecar(wavecar: &Wavecar) -> (Array3<f64>, Array3<f64>, Array2<f64>) {
        let eigs = wavecar.band_eigs.clone();
        let occs = wavecar.band_fweights.clone();
        let kvec = wavecar.kvecs.clone();

        (eigs, occs, kvec)
    }

    //                                       eigs         occs         kvec
    fn get_bands_kpoints(&self) -> Result<(Array3<f64>, Array3<f64>, Array2<f64>)> {
        Wavecar::from_file(&self.wavecar).and_then(|v| {
            info!("Trying to parse {:?} ...", self.wavecar);
            Ok(Self::bands_from_wavecar(&v))
        })
            .or_else(|_| {
                Procar::from_file(&self.procar).and_then(|v| {
                    info!("Trying to parse {:?} ...", self.procar);
                    Ok(Self::bands_from_procar(&v))
                })
            })
            .with_context(|| format!("Neither WAVECAR nor PROCAR is accessible, please specify a valid WAVECAR or PROCAR"))
    }
}


impl OptProcess for Gap {
    fn process(&self) -> Result<()> {
        let (eigs, occs, kvec) = self.get_bands_kpoints()?;
        let nspin = occs.shape()[0];
        let nkpts = occs.shape()[1];

        let threshold: f64 = occs.iter().copied().fold(f64::NAN, f64::max) / 2.0;

        // lowest conduction band indices
        let cbidx: Array2<usize> = occs.lanes(Axis(2))
            .into_iter()
            .map(|v| v.as_slice().unwrap().partition_point(|&x| x > threshold))
            .collect::<Array1<usize>>()
            .into_shape((nspin, nkpts)).unwrap();

        let vbidx = cbidx.clone() - 1;

        // check if all the cband index are consistent
        let cbi = cbidx[(0, 0)];
        if cbidx.iter().copied().any(|x| x != cbi) {
            let mut output = String::with_capacity(60);
            output.push_str("----------------------------------------\n");
            output.push_str(&format!(" Current system is  {:^20}\n", "Metal".bright_yellow()));
            output.push_str("----------------------------------------");
            println!("{}", output);

            return Ok(());
        }

        // find cbm
        let cbeigs = multizip((cbidx.clone(), eigs.lanes(Axis(2))))
            .map(|(i, v)| v[i])
            .collect::<Array1<f64>>()
            .into_shape((nspin, nkpts)).unwrap();
        let cbm = cbeigs.rows().into_iter()
            .map(|v| v.iter().copied().fold(f64::NAN, f64::min))
            .collect::<Vec<f64>>();
        let cbmik = cbeigs.rows().into_iter().zip(cbm.iter())
            .map(|(v, c)| v.iter().position(|x| x == c).unwrap())
            .collect::<Vec<usize>>();


        // find vbm
        let vbeigs = multizip((vbidx, eigs.lanes(Axis(2))))
            .map(|(i, v)| v[i])
            .collect::<Array1<f64>>()
            .into_shape((nspin, nkpts)).unwrap();
        let vbm = vbeigs.rows()
            .into_iter()
            .map(|v| v.iter().copied().fold(f64::NAN, f64::max))
            .collect::<Vec<f64>>();
        let vbmik = vbeigs.rows().into_iter().zip(vbm.iter())
            .map(|(v, c)| v.iter().position(|x| x == c).unwrap())
            .collect::<Vec<usize>>();

        
        let mut output = String::with_capacity(60);
        output.push_str("--------------------------------------------------------------------------------\n");
        if 1 == nspin {
            let is_direct = if cbmik[0] == vbmik[0] { "Direct Gap" } else { "Indirect Gap" };
            let gap = cbm[0] - vbm[0];
            let kcbm = kvec.row(cbmik[0]);
            let kvbm = kvec.row(vbmik[0]);

            output.push_str(&format!(" Current system has {:^16} of {:^10} eV\n",
                                     is_direct.bright_yellow(), format!("{:5.3}", gap).bright_cyan()
                                     ));
            output.push_str(&format!("  CBM @ k-point {:5} of ({:6.3},{:6.3},{:6.3}) , band {:5} of {:8} eV\n",
                                     cbmik[0], kcbm[0], kcbm[1], kcbm[2], cbidx[(0,0)],
                                     format!("{:8.3}", cbeigs[(0, cbmik[0])]).bright_blue()
                                     ));
            output.push_str(&format!("  VBM @ k-point {:5} of ({:6.3},{:6.3},{:6.3}) , band {:5} of {:8} eV\n",
                                     vbmik[0], kvbm[0], kvbm[1], kvbm[2], cbidx[(0,0)],
                                     format!("{:8.3}", vbeigs[(0, vbmik[0])]).bright_blue()
                                     ));
        } else {
            let spin_ud = ["SPIN UP", "SPIN DOWN"];
            for ispin in 0 .. nspin {
                output.push_str(&format!("    ====================  Gap Info For {:^15}  ====================    \n", spin_ud[ispin].bright_green()));

                let is_direct = if cbmik[ispin] == vbmik[ispin] { "Direct Gap" } else { "Indirect Gap" };
                let gap = cbm[ispin] - vbm[ispin];
                let kcbm = kvec.row(cbmik[ispin]);
                let kvbm = kvec.row(vbmik[ispin]);

                output.push_str(&format!(" Current system has {:^16} of {:^10} eV\n",
                                         is_direct.bright_yellow(), format!("{:5.3}", gap).bright_cyan()
                                         ));
                output.push_str(&format!("  CBM @ k-point {:5} of ({:6.3},{:6.3},{:6.3}) , band {:5} of {:8} eV\n",
                                         cbmik[ispin], kcbm[0], kcbm[1], kcbm[2], cbidx[(ispin,0)],
                                         format!("{:8.3}", cbeigs[(ispin, cbmik[ispin])]).bright_blue()
                                         ));
                output.push_str(&format!("  VBM @ k-point {:5} of ({:6.3},{:6.3},{:6.3}) , band {:5} of {:8} eV\n",
                                         vbmik[ispin], kvbm[0], kvbm[1], kvbm[2], cbidx[(ispin,0)],
                                         format!("{:8.3}", vbeigs[(ispin, vbmik[ispin])]).bright_blue()
                                         ));
            }
        }
        output.push_str("--------------------------------------------------------------------------------");


        println!("{}", output);

        Ok(())
    }
}
