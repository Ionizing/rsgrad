use std::{
    path::PathBuf,
    fs,
};
use log::{
    info,
    warn,
    debug,
};
use clap::Args;
use crate::{
    Result,
    OptProcess,
    Outcar,
    IonicIterationsFormat,
    Poscar,
};


#[derive(Debug, Args)]
/// Tracking relaxation or MD progress.
///
/// Contains the evolution of energy, maximum of Hellmann-Feynman forces, 
/// magnetic moments and time usage of each ionic step.
///
/// Hint: This command may require POSCAR for atom constraints information.
pub struct Rlx {
    #[arg(default_value = "./OUTCAR")]
    /// Specify the input OUTCAR file
    outcar: PathBuf,

    #[arg(short = 'p', long, default_value = "./POSCAR")]
    /// Specify the input POSCAR file
    poscar: PathBuf,

    #[arg(short = 'e', long = "toten")]
    /// Prints TOTEN in eV
    print_energy: bool,

    #[arg(short = 'a', long = "favg")]
    /// Prints averaged total force in eV/A
    print_favg: bool,

    #[arg(short = 'x', long = "fmaxis")]
    /// Prints the axis where the strongest total force component lies on. [XYZ]
    print_fmax_axis: bool,

    #[arg(short = 'i' ,long = "fmidx")]
    /// Prints the index of ion with maximum total force load. Starts from 1
    print_fmax_index: bool,

    #[arg(short = 'v', long = "volume")]
    /// Prints lattice volume in A^3
    print_volume: bool,

    #[arg(long = "no-fmax")]
    /// Don't print maximum total force in A^3
    no_print_fmax: bool,

    #[arg(long = "no-totenz")]
    /// Don't print TOTEN without entropy in eV
    no_print_energyz: bool,

    #[arg(long = "no-lgde")]
    /// Don't print Log10(delta(TOTEN without entropy))
    no_print_lgde: bool,

    #[arg(long = "no-magmom")]
    /// Don't print total magnetic moment in muB
    no_print_magmom: bool,

    #[arg(long = "no-nscf")]
    /// Don't print number of SCF iteration for each ionic step
    no_print_nscf: bool,

    #[arg(long = "no-time")]
    /// Don't print time elapsed for each ionic step in minutes
    no_print_time: bool,
}


impl OptProcess for Rlx {
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

        let iif = IonicIterationsFormat::from_outcar(&outcar)
            .print_energy     ( self.print_energy)
            .print_energyz    (!self.no_print_energyz)
            .print_log10de    (!self.no_print_lgde)
            .print_favg       ( self.print_favg)
            .print_fmax       (!self.no_print_fmax)
            .print_fmax_axis  ( self.print_fmax_axis)
            .print_fmax_index ( self.print_fmax_index)
            .print_nscf       (!self.no_print_nscf)
            .print_time_usage (!self.no_print_time)
            .print_magmom     (!self.no_print_magmom)
            .print_volume     ( self.print_volume);
        print!("{}", iif);
        Ok(())
    }
}
