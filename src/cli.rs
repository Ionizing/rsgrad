use clap::{
    AppSettings,
    Parser,
};
use enum_dispatch::enum_dispatch;

use crate::{
    types::Result,
    commands::{
        rlx::Rlx,
        vib::Vib,
        traj::Traj,
        pos::Pos,
        pot::Pot,
        chgdiff::Chgdiff,
        workfunc::Workfunc,
        dos::Dos,
        band::Band,
        wav3d::Wav3D,
        wav1d::Wav1D,
        tdm::Tdm,
        gap::Gap,
    },
};

#[enum_dispatch]
pub trait OptProcess : Parser {
    fn process(&self) -> Result<()>;
}

#[enum_dispatch(OptProcess)]
#[derive(Debug, Parser)]
#[clap(name = "rsgrad",
            about = r"A command-line tool to help VASP players play better with VASP.
If you want more detailed documentation, just visit https://ionizing.github.io/rsgrad", 
            version,
            author = "@Ionizing github.com/Ionizing/rsgrad",
            setting = AppSettings::ColoredHelp,
            setting = AppSettings::ColorAuto,
            )]
enum Opt {
    #[clap(setting = AppSettings::ColoredHelp,
           setting = AppSettings::ColorAuto)]
    Rlx,


    #[clap(setting = AppSettings::ColoredHelp,
           setting = AppSettings::ColorAuto)]
    Vib,


    #[clap(setting = AppSettings::ColoredHelp,
           setting = AppSettings::ColorAuto)]
    Traj,


    #[clap(setting = AppSettings::ColoredHelp,
           setting = AppSettings::ColorAuto)]
    Pos,


    #[clap(setting = AppSettings::ColoredHelp,
           setting = AppSettings::ColorAuto)]
    Pot,


    #[clap(setting = AppSettings::ColoredHelp,
           setting = AppSettings::ColorAuto)]
    Chgdiff,


    #[clap(setting = AppSettings::ColoredHelp,
           setting = AppSettings::ColorAuto)]
    Workfunc,


    #[clap(setting = AppSettings::ColoredHelp,
           setting = AppSettings::ColorAuto)]
    Dos,


    #[clap(setting = AppSettings::ColoredHelp,
           setting = AppSettings::ColorAuto)]
    Band,


    #[clap(setting = AppSettings::ColoredHelp,
           setting = AppSettings::ColorAuto,
           name = "wav3d")]
    Wav3D,


    #[clap(setting = AppSettings::ColoredHelp,
           setting = AppSettings::ColorAuto,
           name = "wav1d")]
    Wav1D,


    #[clap(setting = AppSettings::ColoredHelp,
           setting = AppSettings::ColorAuto)]
    Tdm,


    #[clap(setting = AppSettings::ColoredHelp,
           setting = AppSettings::ColorAuto)]
    Gap
}


pub fn run() -> Result<()> {
    Opt::parse().process()
}
