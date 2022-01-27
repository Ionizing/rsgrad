use anyhow::{
    Context,
    Result,
};
use vasp_poscar::Poscar;


pub struct AtomPotcar {
    pub symbol: String,  // Element symbol, H, He, Li, Be, B, C ...
    pub functional: String,  // functional type, GGA, LDA, PBE, GW ...
    pub content: String,  // raw content of single element POTCAR
}


impl AtomPotcar {
    fn from_config() -> Self {
        todo!();
    }
}


pub struct Potcar {
    pub data: Vec<AtomPotcar>,
}
