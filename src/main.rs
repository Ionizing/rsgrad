use std::io::Result;
use std::path::PathBuf;
use std::time;
use env_logger;
use log::{
    info,
    warn,
    debug,
};
use structopt::StructOpt;
use rsgrad::outcar::Outcar;
use rsgrad::format::{
    IonicIterationsFormat,
    Vibrations,
    PrintAllVibFreqs,
};

use structopt::clap::AppSettings;

#[derive(Debug, StructOpt)]
#[structopt(name = "rsgrad",
            about = "A tool used to tracking the relaxation or MD progress of VASP calculation",
            author = "@Ionizing  <PeterSmith_9@outlook.com>",
            setting = AppSettings::ColoredHelp,
            setting = AppSettings::ColorAuto)]
struct Opt {
    #[structopt(subcommand)]
    command: Command,

    #[structopt(default_value = "./OUTCAR")]
    /// Specify the input OUTCAR file name
    input: PathBuf,
}

#[derive(Debug, StructOpt)]
enum Command {
    #[structopt(setting = AppSettings::ColoredHelp,
                setting = AppSettings::ColorAuto)]
    /// Tracking info associated with relaxation stuff
    Rlx {
        #[structopt(short = "e", long = "toten")]
        /// Prints TOTEN in eV
        print_energy: bool,

        #[structopt(short = "a", long = "favg")]
        /// Prints averaged total force in eV/A
        print_favg: bool,

        #[structopt(short = "x", long = "fmaxis")]
        /// Prints the axis where the strongest total force component lies on. [XYZ]
        print_fmax_axis: bool,

        #[structopt(short = "i" ,long = "fmidx")]
        /// Prints the index of ion with maximum total force load. Starts from 1
        print_fmax_index: bool,

        #[structopt(short = "v", long = "volume")]
        /// Prints lattice volume in A^3
        print_volume: bool,

        #[structopt(long = "no-fmax")]
        /// Don't print maximum total force in A^3
        no_print_fmax: bool,

        #[structopt(long = "no-totenz")]
        /// Don't print TOTEN without entropy in eV
        no_print_energyz: bool,

        #[structopt(long = "no-lgde")]
        /// Don't print Log10(delta(TOTEN without entropy))
        no_print_lgde: bool,

        #[structopt(long = "no-magmom")]
        /// Don't print total magnetic moment in muB
        no_print_magmom: bool,

        #[structopt(long = "no-nscf")]
        /// Don't print number of SCF iteration for each ionic step
        no_print_nscf: bool,

        #[structopt(long = "no-time")]
        /// Don't print time elapsed for each ionic step in minutes
        no_print_time: bool,
    },

    #[structopt(setting = AppSettings::ColoredHelp,
                setting = AppSettings::ColorAuto,
                setting = AppSettings::AllowNegativeNumbers)]
    /// Tracking info associated with vibration stuff
    Vib {
        #[structopt(short, long)]
        /// Shows vibration modes in brief
        list: bool,

        #[structopt(short = "x", long)]
        save_as_xsfs: bool,

        #[structopt(short = "i", long)]
        select_indices: Option<Vec<i32>>,

        #[structopt(long, default_value = ".")]
        save_in: PathBuf,
    },

    #[structopt(setting = AppSettings::ColoredHelp,
                setting = AppSettings::ColorAuto,
                setting = AppSettings::AllowNegativeNumbers)]
    /// Operations about relaxation/MD trajectory
    Trj {
        #[structopt(short, long)]
        list: bool,

        #[structopt(short = "i", long)]
        /// Select the indices to operate.
        ///
        /// Step indices start from '1', if '0' is given, all the structures will be selected.
        /// Step indices can be negative, where negative index means counting reversely.
        /// E.g. "--save-as-poscars -2 -1 1 2 3" means saving the last two and first three
        /// steps.
        select_indices: Option<Vec<i32>>,

        #[structopt(short = "d", long)]
        /// Save total trajectory in XDATCAR format
        save_as_xdatcar: bool,

        #[structopt(short = "p", long)]
        /// Save structures of given steps as POSCARs.
        save_as_poscars: bool,

        #[structopt(short = "x", long)]
        /// Save structures of given steps as XSFs.
        save_as_xsfs: bool,

        #[structopt(long, default_value = ".")]
        /// Define where the files would be saved.
        save_in: PathBuf,
    },

}


fn main() -> Result<()> {
    let now = time::Instant::now();

    let env = env_logger::Env::new().filter_or("RSGRAD_LOG", "info");
    env_logger::init_from_env(env);

    let opt = Opt::from_args();
    debug!("{:?}", opt);

    info!("Parsing input file {:?} ...", &opt.input);
    let f = Outcar::from_file(&opt.input)?;

    match opt.command {
        Command::Rlx { print_energy,
                       print_favg,
                       print_fmax_axis,
                       print_fmax_index,
                       print_volume,
                       no_print_fmax,
                       no_print_energyz,
                       no_print_lgde,
                       no_print_magmom,
                       no_print_nscf,
                       no_print_time } => {
            let iif = IonicIterationsFormat::from(f.ion_iters)
                .print_energy     (print_energy)
                .print_energyz    (!no_print_energyz)
                .print_log10de    (!no_print_lgde)
                .print_favg       (print_favg)
                .print_fmax       (!no_print_fmax)
                .print_fmax_axis  (print_fmax_axis)
                .print_fmax_index (print_fmax_index)
                .print_nscf       (!no_print_nscf)
                .print_time_usage (!no_print_time)
                .print_magmom     (!no_print_magmom)
                .print_volume     (print_volume);
            print!("{}", iif);
        },
        Command::Vib { list,
                       save_as_xsfs,
                       select_indices,
                       save_in } => {
            if list {
                let paf: PrintAllVibFreqs = Vibrations::from(f).into();
                print!("{}", paf);
                return Ok(());
            }

            if save_as_xsfs {
                let select_indices = select_indices.unwrap_or_default();
                if select_indices.len() == 0 {
                    warn!("No modes are selected to operate!");
                    return Ok(());
                }

                let vibs = Vibrations::from(f);
                let len = vibs.modes.len();

                let inds: Vec<usize> =
                    if select_indices.contains(&0) {
                        (1..=len).collect()
                    } else {
                        select_indices
                            .into_iter()
                            .map(|i| {
                                if i < 0 {
                                    i.rem_euclid(len as i32) as usize + 1
                                } else {
                                    i as usize + 1
                                }
                            })
                            .collect::<Vec<usize>>()
                    };

                for i in inds {
                    info!("Saving mode #{:4} into {:?} ...", i, save_in);
                    vibs.save_as_xsf(i, &save_in)?;
                }
            }
        },
        Command::Trj { .. } => {

        }
    }

    info!("Time used: {:?}", now.elapsed());
    Ok(())
}
