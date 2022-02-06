use std::time;

use env_logger;
use log::info;
use rsgrad::traits::Result;
use structopt::StructOpt;
use structopt::clap::AppSettings;
use rsgrad::traits::OptProcess;
use rsgrad::commands;

#[derive(Debug, StructOpt)]
#[structopt(name = "rsgrad",
            about = "A tool used to track the VASP calculation result",
            author = "@Ionizing github.com/Ionizing/rsgrad",
            setting = AppSettings::ColoredHelp,
            setting = AppSettings::ColorAuto,
            )]
enum Opt {
    #[structopt(setting = AppSettings::ColoredHelp,
                setting = AppSettings::ColorAuto)]
    Rlx(commands::rlx::Rlx),


    #[structopt(setting = AppSettings::ColoredHelp,
                setting = AppSettings::ColorAuto)]
    Vib(commands::vib::Vib),


    #[structopt(setting = AppSettings::ColoredHelp,
                setting = AppSettings::ColorAuto)]
    Traj(commands::traj::Traj),


    #[structopt(setting = AppSettings::ColoredHelp,
                setting = AppSettings::ColorAuto)]
    Pos(commands::pos::Pos),
}


fn main() -> Result<()> {
    let now = time::Instant::now();

    env_logger::init_from_env(
        env_logger::Env::new().filter_or("RSGRAD_LOG", "info"));

    match Opt::from_args() {
        Opt::Rlx(cmd) => cmd.process()?,
        Opt::Vib(cmd) => cmd.process()?,
        Opt::Traj(cmd) => cmd.process()?,
        Opt::Pos(cmd) => cmd.process()?,
    }

    info!("Time used: {:?}", now.elapsed());
    Ok(())
}
