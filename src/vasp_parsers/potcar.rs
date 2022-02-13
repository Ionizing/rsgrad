use anyhow::{
    Context,
    Result,
};


pub enum FunctionalType {
    PawPbe,
    PawLda,
}


pub struct AtomicPotcar {
    pub symbol: String,                 // Element symbol, H, He, Li, Be, B, C ...
    pub functional: FunctionalType,     // Functional type, LDA, PBE
    pub valence_type: Vec<String>,      // Valence annotation '_sv', '_GW', '_AE' ...
    pub content: String,                // Raw content of single element POTCAR
}


impl AtomicPotcar {
    fn from_config() -> Self {
        todo!();
    }
}


pub struct Potcar {
    pub data: Vec<AtomicPotcar>,
}
