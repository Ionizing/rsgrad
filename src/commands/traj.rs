use std::{
    path::PathBuf,
    fs,
};
use log::{
    info,
    warn,
    debug,
};
use rayon::prelude::*;
use structopt::{
    StructOpt,
    clap::AppSettings,
};
use crate::{
    Result,
    index_transform,
    OptProcess,
    Outcar,
    Poscar,
    Trajectory,
};


#[derive(Debug, StructOpt)]
#[structopt(setting = AppSettings::ColoredHelp,
            setting = AppSettings::ColorAuto,
            setting = AppSettings::AllowNegativeNumbers)]
/// Operations about relaxation/MD trajectory.
///
/// POSCAR is needed if you want to preserve the constraints when saving frames to POSCAR.
pub struct Traj {
    #[structopt(default_value = "./OUTCAR")]
    /// Specify the input OUTCAR file
    outcar: PathBuf,

    #[structopt(short = "p", long, default_value = "./POSCAR")]
    /// Specify the input POSCAR file
    poscar: PathBuf,

    #[structopt(short = "x", long)]
    /// Saves each selected modes to XSF file, this file includes each atom's force information
    save_as_xsfs: bool,

    #[structopt(short = "s", long)]
    /// Save selected steps as POSCARs
    save_as_poscar: bool,

    #[structopt(short = "d", long)]
    /// Save whole trajectory in XDATCAR format
    save_as_xdatcar: bool,

    #[structopt(short = "i", long)]
    /// Selects the indices to operate.
    ///
    /// Step indices start from '1', if '0' is given, all the structures will be selected.
    /// Step indices can be negative, where negative index means counting reversely.
    /// E.g. "--save-as-poscars -2 -1 1 2 3" means saving the last two and first three
    /// steps.
    select_indices: Option<Vec<i32>>,

    #[structopt(long, default_value = ".")]
    /// Define where the files would be saved
    save_in: PathBuf,

    #[structopt(long = "no-add-symbol-tags")]
    /// Don't add chemical symbol to each line of coordinates
    no_add_symbol_tags: bool,

    #[structopt(long = "no-preserve-constraints")]
    /// Don't preverse constraints when saving trajectory to POSCAR
    no_preserve_constraints: bool,

    #[structopt(long = "cartesian")]
    /// Save to POSCAR in cartesian coordinates, the coordinates written is direct/fractional by
    /// default
    cartesian: bool,
}


impl OptProcess for Traj {
    fn process(&self) -> Result<()> {
        info!("Parsing file {:?} and {:?}", &self.outcar, &self.poscar);
        debug!("    OUTCAR file path = {:?}\n    POSCAR file path = {:?}",
               fs::canonicalize(&self.outcar), fs::canonicalize(&self.poscar));

        let mut outcar = Outcar::from_file(&self.outcar)?;
        if let Ok(poscar) = Poscar::from_file(&self.poscar) {
            if let Some(constraints) = poscar.constraints {
                outcar.set_constraints(constraints);
            }
        } else {
            warn!("Reading constraints from POSCAR file {:?} failed", &self.poscar);
        }

        let traj = Trajectory::from(outcar.clone());

        if self.save_as_xdatcar {
            traj.save_as_xdatcar(&self.save_in)?;
        }

        let inds = {
            let select_indices = self.select_indices.clone().unwrap_or_default();
            if 0 == select_indices.len() {
                warn!("No steps are selected to operate !");
            }
            index_transform(select_indices, traj.0.len())
        };

        if self.save_as_poscar {
            inds.par_iter()
                .map(|i| {
                    traj.save_as_poscar(*i, &self.save_in, 
                                        !self.cartesian, 
                                        !self.no_preserve_constraints, 
                                        !self.no_add_symbol_tags)?;
                    Ok(())
                })
            .collect::<Result<()>>()?;
        }

        if self.save_as_xsfs {
            inds.par_iter()
                .map(|i| {
                    outcar.save_ionic_step_as_xsf(*i, &self.save_in)?;
                    Ok(())
                })
            .collect::<Result<()>>()?;
        }

        Ok(())
    }
}
