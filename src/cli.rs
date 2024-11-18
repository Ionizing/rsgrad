use std::sync::OnceLock;
use clap::{
    Parser,
    builder::styling::{
        AnsiColor,
        Effects,
        Styles,
    },
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
        uc::Uc, 
        modelnac::ModelNac,
    },
};


pub fn get_style() -> Styles {
    static INSTANCE: OnceLock<Styles> = OnceLock::new();
    INSTANCE.get_or_init(|| {
        Styles::styled()
            .header(AnsiColor::Yellow.on_default() | Effects::BOLD)
            .usage(AnsiColor::Green.on_default()   | Effects::BOLD)
            .literal(AnsiColor::Green.on_default() | Effects::BOLD)
            .placeholder(AnsiColor::BrightBlue.on_default())
            .error(AnsiColor::BrightRed.on_default())
            .valid(AnsiColor::BrightYellow.on_default())
    }).to_owned()
}


#[enum_dispatch]
pub trait OptProcess {
    fn process(&self) -> Result<()>;
}


#[enum_dispatch(OptProcess)]
#[derive(Debug, Parser)]
#[command(name = "rsgrad",
            about = r"A command-line tool to help VASP players play better with VASP.
If you want more detailed documentation, just visit https://ionizing.github.io/rsgrad", 
            version,
            author = "@Ionizing github.com/Ionizing/rsgrad",
            styles = get_style()
            )]
enum Opt {
    Rlx,

    Vib,

    Traj,

    Pos,

    Pot,

    Chgdiff,

    Workfunc,

    Dos,
    
    Band,

    #[command(name = "wav3d")]
    Wav3D,

    #[command(name = "wav1d")]
    Wav1D,

    Tdm,

    Gap,

    Uc,

    ModelNac,
}


pub fn run() -> Result<()> {
    Opt::parse().process()
}
