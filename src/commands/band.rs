use std::{
    fs,
    path::PathBuf,
};

use indexmap::IndexMap;
use structopt::{
    StructOpt,
    clap::AppSettings,
};
use log::{
    info,
    warn,
};
use rayon;
use anyhow::{
    bail,
    anyhow,
    Context,
};
use serde::{
    Serialize,
    Deserialize,
};
use ndarray::{
    s,
    arr2,
    Axis,
};
use itertools::Itertools;

use crate::{
    Result,
    OptProcess,
    Procar,
    Outcar,
    Poscar,
    vasp_parsers::outcar::GetEFermi,
    types::{
        Mat33,
        MatX3,
        Vector,
        Matrix,
        Cube,
    },
    commands::common::{
        write_array_to_txt,
        RawSelection,
    }
};


const THRESHOLD: f64 = 1E-6;


#[derive(Clone, Debug)]
struct Selection {
    label:      String,
    ispins:     Vec<usize>,
    ikpoints:   Vec<usize>,
    iatoms:     Vec<usize>,
    iorbits:    Vec<usize>,
    factor:     f64,
}


fn rawsel_to_sel(r: IndexMap<String, RawSelection>, 
                 nspin: usize,
                 is_ncl: bool,
                 nlm: &[String], 
                 nions: usize, 
                 nkpoints: usize) -> Result<Vec<Selection>> {

    let mut sel_vec = vec![];

    for (label, val) in r.into_iter() {
        let ispins      = RawSelection::parse_ispins(   val.spins.as_deref(),   nspin, is_ncl)?;
        let ikpoints    = RawSelection::parse_ikpoints( val.kpoints.as_deref(), nkpoints)?;
        let iatoms      = RawSelection::parse_iatoms(   val.atoms.as_deref(),   nions)?;
        let iorbits     = RawSelection::parse_iorbits(  val.orbits.as_deref(),  nlm)?;
        let factor      = val.factor.unwrap_or(1.0);

        if factor < 0.0 { bail!("The factor cannot be negative."); }

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


#[derive(Clone, Serialize, Deserialize, Debug)]
struct Configuration {
    kpoint_labels: Option<Vec<String>>,

    #[serde(default = "Configuration::procar_default")]
    procar: PathBuf,

    #[serde(default = "Configuration::outcar_default")]
    outcar: PathBuf,

    #[serde(default = "Configuration::txtout_default")]
    txtout: PathBuf,

    #[serde(default = "Configuration::htmlout_default")]
    htmlout: PathBuf,

    //#[serde(default = "Configuration::is_hse_default")]
    //is_hse: bool,

    #[serde(default = "Configuration::segment_ranges_default")]
    segment_ranges: Option<Vec<(usize, usize)>>,

    pband: Option<IndexMap<String, RawSelection>>,
}

impl Configuration {
    pub fn procar_default() -> PathBuf { PathBuf::from("./PROCAR") }
    pub fn outcar_default() -> PathBuf { PathBuf::from("./OUTCAR") }
    pub fn txtout_default() -> PathBuf { PathBuf::from("./band_raw.txt") }
    pub fn htmlout_default() -> PathBuf { PathBuf::from("./band.html") }
    //pub fn is_hse_default() -> bool { false }
    pub fn segment_ranges_default() -> Option<Vec<(usize, usize)>> { None }
}


#[derive(Debug, StructOpt, Clone)]
#[structopt(setting = AppSettings::ColoredHelp,
            setting = AppSettings::ColorAuto)]
pub struct Band {
    #[structopt(short, long)]
    /// Band structure plot configuration file path.
    ///
    /// If left empty, only bare band is calculated. The configuration template
    /// can be generated by `--gen-template` and then you can follow it.
    config: Option<PathBuf>,

    #[structopt(long)]
    /// Generate band structure plot configuration template.
    gen_template: bool,

    #[structopt(short, long)]
    /// Symbols for high symmetry points on the kpoint path.
    kpoint_labels: Option<Vec<String>>,

    //#[structopt(long)]
    // /// Specify the system is calculated via HSE method or not.
    //is_hse: bool,

    #[structopt(long, default_value = "./PROCAR")]
    /// PROCAR path.
    ///
    /// The band level and projected band info are extracted from PROCAR.
    procar: PathBuf,

    #[structopt(long, default_value = "./OUTCAR")]
    /// OUTCAR path.
    ///
    /// The fermi level and lattice info are extracted from OUTCAR.
    outcar: PathBuf,

    #[structopt(long, default_value = "band.txt")]
    /// Save the raw data of band structure.
    ///
    /// Then you can replot it with more advanced tools.
    txtout: PathBuf,

    #[structopt(long, default_value = "band.html")]
    /// Save the band structure plot as HTML.
    ///
    /// Note: Your browser should be able to run plotly.js. Chrome, Safari, Edge, Firefox and
    /// etc. are supported.
    htmlout: PathBuf,
}

// Extra TODO: shift eigvals to E-fermi
impl Band {
    /// Given a `nkpoints*3` matrix to generate k-point path for bandstructure
    /// `segments_ranges` are the start and end indices of each, closed interval.
    /// the range can be in reversed direction. Index starts from 1.
    fn gen_kpath(kpoints: &Matrix<f64>, bcell: &[[f64; 3]; 3], segment_ranges: &Vec<(usize, usize)>)-> Vector<f64> {
        let bcell = arr2(bcell);
        let kpoints = kpoints.dot(&bcell);

        segment_ranges
            .iter()
            .cloned()
            .map(|(beg, end)| {
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
            .flatten()
            .scan(0.0, |acc, x| {
                *acc += x;
                Some(*acc)
            })
            .collect()
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
            if (kpoints.index_axis(Axis(0), i).to_owned() - 
                kpoints.index_axis(Axis(0), i-1))
                .iter().all(|d| d.abs() < THRESHOLD) {  // if kpt[i, :] == kpt[i-1, :]

                boundaries.push((last_bound, i));
                last_bound = i + 1;
            }
        }

        boundaries.push((last_bound, len));
        Ok(boundaries)
    }


    /// Return: [ispin, ikpoint, iband]
    fn plot_rawband(eigvals: Cube<f64>) -> Vec<Vec<f64>> {
        todo!()
    }

    fn plot_pband() ->Vec<f64> {
        todo!()
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


impl OptProcess for Band {
    fn process(&self) -> Result<()> {
        let mut procar: Result<Procar> = Err(anyhow!(""));
        let mut outcar: Result<Outcar> = Err(anyhow!(""));

        rayon::scope(|s| {
            s.spawn(|_| {
                info!("Reading band data from {:?}", &self.procar);
                procar = Procar::from_file(&self.procar);
            });
            s.spawn(|_| {
                info!("Reading fermi level and lattice data from {:?}", &self.outcar);
                outcar = Outcar::from_file(&self.outcar);
            });
        });

        let mut procar = procar.context(format!("Parse PROCAR file {:?} failed.", self.procar))?;
        let outcar = outcar.context(format!("Parse OUTCAR file {:?} failed.", self.procar))?;

        let efermi = outcar.efermi;
        let cell = outcar.ion_iters.last()
            .context("This OUTCAR doesn't complete at least one ionic step.")?
            .cell;

        info!("Found Fermi level: {}, shifting eigenvalues ...", efermi);
        procar.pdos.eigvals -= efermi;

        //if self.is_hse {
            //if Self::filter_hse(&mut procar) {
                //info!("Zero-contribution k-points filtered out for HSE system.")
            //} else {
                //warn!("Could not find zero-contribution k-points, processing the non-zero-contribution k-points.")
            //}
        //}
        
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

        let cell = [[    2.8847306200059699 * 0.993,  -1.6655000000000000 * 0.993,   0.0000000000000000 * 0.993],
                    [    0.0000000000000000 * 0.993,   3.3310000000000000 * 0.993,   0.0000000000000000 * 0.993],
                    [    0.0000000000000000 * 0.993,   0.0000000000000000 * 0.993,  23.1640000000000015 * 0.993]];

        let bcell = Poscar::acell_to_bcell(&cell).unwrap();

        let kpath = Band::gen_kpath(&kpoints, &bcell, &vec![(1, 5), (6, 10), (11, 15), (6, 10)]);
        let expect = arr1(&[0.0, 0.05041295144327735, 0.1008259028865547, 0.15123885432983203, 0.2016518057731094,
                          0.2016518057731094, 0.22682051907635853, 0.2519892323796077, 0.27715794568285684, 0.30232665898610595,
                          0.30232665898610595, 0.3459637207291467, 0.38960078247218755, 0.4332378442152283, 0.4768749059582691,
                          0.4768749059582691, 0.5020436192615182, 0.5272123325647673, 0.5523810458680164, 0.5775497591712655]);
        eprintln!("{}", &kpath);
        assert!((kpath - &expect).iter().all(|x| x.abs() < 1E-6));

        let segment_ranges = Band::find_segments(&kpoints).unwrap();
        assert_eq!(segment_ranges, vec![(1, 5), (6, 10), (11, 15)]);
    }
}
