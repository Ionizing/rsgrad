use std::time;

use env_logger;
use log::info;
use rsgrad::{
    Result,
    OptProcess,
    commands,
};
use structopt::{
    clap::AppSettings,
    StructOpt,
};

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


    #[structopt(setting = AppSettings::ColoredHelp,
                setting = AppSettings::ColorAuto)]
    Pot(commands::pot::Pot),


    #[structopt(setting = AppSettings::ColoredHelp,
                setting = AppSettings::ColorAuto)]
    Chgdiff(commands::chgdiff::Chgdiff),


    #[structopt(setting = AppSettings::ColoredHelp,
                setting = AppSettings::ColorAuto)]
    Workfunc(commands::workfunc::Workfunc),
}


fn main() -> Result<()> {
    let now = time::Instant::now();

    env_logger::init_from_env(
        env_logger::Env::new().filter_or("RSGRAD_LOG", "info"));

    match Opt::from_args() {
        Opt::Rlx(cmd)       => cmd.process()?,
        Opt::Vib(cmd)       => cmd.process()?,
        Opt::Traj(cmd)      => cmd.process()?,
        Opt::Pos(cmd)       => cmd.process()?,
        Opt::Pot(cmd)       => cmd.process()?,
        Opt::Chgdiff(cmd)   => cmd.process()?,
        Opt::Workfunc(cmd)  => cmd.process()?,
    }

    info!("Time used: {:?}", now.elapsed());
    Ok(())
}
