use std::{
    fmt,
    fs::File
};

use crate::types::Axis;

/// Wavefunction precision type
///
/// Usually VASP stores the band coefficients in two precisions: float32 and float64.
/// From the precision tag in WAVECAR's header we can infer the precision type.
#[derive(PartialEq, Debug, Clone, Copy)]
pub enum WFPrecType {
    Complex32,
    Complex64,
}

/// WAVECAR type enumeration
#[derive(Debug, Clone, Copy)]
pub enum WavecarType {
    /// Most typical WAVECAR type. VASP is executed via `vasp_std`, calculation is done
    /// at all k-points.
    Standard,
    /// VASP is executed via `vasp_gam`, only half of the k-space is utilized because
    /// only Gamma point is sampled, which meas the wavefunction in real space is totally
    /// real, thus only half of the coefficients are needed in k-space.
    GamaHalf(Axis),
    /// VASP is executed via `vasp_ncl`, there should be `LNONCOLLINEAR = T` in OUTCAR.
    /// In this type of WAVECAR, ISPIN = 1, but two spinor components are stored in ONE
    /// spin component side by side.
    NonCollinear,
}

impl fmt::Display for WavecarType {
    fn fmt(&self, f: &mut fmt::Formatter<'_>) -> fmt::Result {
        let description = match self {
            Self::Standard          => "Standard",
            Self::NonCollinear      => "NonCollinear",
            Self::GamaHalf(Axis::X) => "GammaX",
            Self::GamaHalf(Axis::Y) => "GammaY",
            Self::GamaHalf(Axis::Z) => "GammaZ",
        };
        f.write_str(description)
    }
}


pub struct Wavecar {

}
