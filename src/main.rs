use std::io::Result;
use std::path::PathBuf;
use std::time;
use log::{
    info,
    Record,
    LevelFilter,
    Level,
    Metadata
};
use structopt::StructOpt;
use rsgrad::outcar::Outcar;
use rsgrad::format::IonicIterationsFormat;

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

    #[structopt(short, long)]
    /// Show total time usage
    time: bool,
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
                setting = AppSettings::ColorAuto)]
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


struct Logger;
static LOGGER: Logger = Logger;

impl log::Log for Logger {
    fn enabled(&self, metadata: &Metadata) -> bool {
        metadata.level() <= Level::Info
    }

    fn log(&self, record: &Record) {
        if self.enabled(record.metadata()) {
            println!("# [{}] -- {}", record.level(), record.args());
        }
    }
    fn flush(&self) { }
}


fn main() -> Result<()> {
    let now = time::Instant::now();

    // set logger
    log::set_logger(&LOGGER).unwrap();
    log::set_max_level(LevelFilter::Info);

    let opt = Opt::from_args();

    info!("{:?}", opt);

    // let f = Outcar::from_file(Path::new(matches.value_of("input").unwrap()))?;

    // let iif = IonicIterationsFormat::from(f.ion_iters)
    //     .print_energy     (convert_helper(matches.value_of("print-energy").unwrap()))
    //     .print_energyz    (convert_helper(matches.value_of("print-energyz").unwrap()))
    //     .print_log10de    (convert_helper(matches.value_of("print-lgde").unwrap()))
    //     .print_favg       (convert_helper(matches.value_of("print-favg").unwrap()))
    //     .print_fmax       (convert_helper(matches.value_of("print-fmax").unwrap()))
    //     .print_fmax_axis  (convert_helper(matches.value_of("print-fmax-axis").unwrap()))
    //     .print_fmax_index (convert_helper(matches.value_of("print-fmax-index").unwrap()))
    //     .print_nscf       (convert_helper(matches.value_of("print-nscf").unwrap()))
    //     .print_time_usage (convert_helper(matches.value_of("print-time").unwrap()))
    //     .print_magmom     (convert_helper(matches.value_of("print-magmom").unwrap()))
    //     .print_volume     (convert_helper(matches.value_of("print-volume").unwrap()));

    // print!("{}", iif);

    info!("Time used: {:?}", now.elapsed());
    Ok(())
}
