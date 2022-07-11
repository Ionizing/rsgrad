use std::time;

use env_logger::init_from_env;
use log::info;
use rsgrad::{
    Result,
    OptProcess,
    commands,
};
use clap::{
    Parser,
    AppSettings,
};

#[derive(Debug, Parser)]
#[clap(name = "rsgrad",
            about = "A tool used to track the VASP calculation result",
            author = "@Ionizing github.com/Ionizing/rsgrad",
            setting = AppSettings::ColoredHelp,
            setting = AppSettings::ColorAuto,
            )]
enum Opt {
    #[clap(setting = AppSettings::ColoredHelp,
                setting = AppSettings::ColorAuto)]
    Rlx(commands::rlx::Rlx),


    #[clap(setting = AppSettings::ColoredHelp,
                setting = AppSettings::ColorAuto)]
    Vib(commands::vib::Vib),


    #[clap(setting = AppSettings::ColoredHelp,
                setting = AppSettings::ColorAuto)]
    Traj(commands::traj::Traj),


    #[clap(setting = AppSettings::ColoredHelp,
                setting = AppSettings::ColorAuto)]
    Pos(commands::pos::Pos),


    #[clap(setting = AppSettings::ColoredHelp,
                setting = AppSettings::ColorAuto)]
    Pot(commands::pot::Pot),


    #[clap(setting = AppSettings::ColoredHelp,
                setting = AppSettings::ColorAuto)]
    Chgdiff(commands::chgdiff::Chgdiff),


    #[clap(setting = AppSettings::ColoredHelp,
                setting = AppSettings::ColorAuto)]
    Workfunc(commands::workfunc::Workfunc),


    #[clap(setting = AppSettings::ColoredHelp,
                setting = AppSettings::ColorAuto)]
    Dos(commands::dos::Dos),


    #[clap(setting = AppSettings::ColoredHelp,
                setting = AppSettings::ColorAuto)]
    Band(commands::band::Band),


    #[clap(setting = AppSettings::ColoredHelp,
                setting = AppSettings::ColorAuto,
                name = "wav3d")]
    Wav3D(commands::wav3d::Wav3D),


    #[clap(setting = AppSettings::ColoredHelp,
                setting = AppSettings::ColorAuto,
                name = "wav1d")]
    Wav1D(commands::wav1d::Wav1D),
}


fn main() -> Result<()> {
    let now = time::Instant::now();

    init_from_env(
        env_logger::Env::new().filter_or("RSGRAD_LOG", "info"));

    match Opt::from_args() {
        Opt::Rlx(cmd)       => cmd.process()?,
        Opt::Vib(cmd)       => cmd.process()?,
        Opt::Traj(cmd)      => cmd.process()?,
        Opt::Pos(cmd)       => cmd.process()?,
        Opt::Pot(cmd)       => cmd.process()?,
        Opt::Chgdiff(cmd)   => cmd.process()?,
        Opt::Workfunc(cmd)  => cmd.process()?,
        Opt::Dos(cmd)       => cmd.process()?,
        Opt::Band(cmd)      => cmd.process()?,
        Opt::Wav3D(cmd)     => cmd.process()?,
        Opt::Wav1D(cmd)     => cmd.process()?,
    }

    info!("Time used: {:?}", now.elapsed());
    Ok(())
}
