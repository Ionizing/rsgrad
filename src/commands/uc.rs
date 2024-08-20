/// Metrix prefixes
pub enum MetricPrefix {
    /// 10⁻¹⁸
    Atto,

    /// 10⁻¹⁵
    Femto,

    /// 10⁻¹²
    Pico,

    /// 10⁻⁹
    Nano,

    /// 10⁻⁶
    Micro,

    /// 10⁻³
    Milli,

    /// 1
    Unity,

    /// 10⁺³
    Kilo,

    /// 10⁺⁶
    Mega,

    /// 10⁺⁹
    Giga,

    /// 10⁺¹²
    Tera,

    /// 10⁺¹⁵
    Peta,

    /// 10⁺¹⁸
    Exa,
}


/// Energy units.
pub enum Unit {
    /// eV, treated as the basic unit
    ElectronVolt,

    /// Cal·mol⁻¹
    CaloriePerMole,

    /// J·mol⁻¹
    JoulePerMole,

    /// Temperature as energy via E=kB*T
    Kelvin,

    /// 1 Hartree ~= 27.2114 eV
    Hartree,

    /// Inverse of wavelength in cm, 1 eV ~= 8065 cm⁻¹
    Wavenumber,

    /// Wavelength of light
    Meter,

    /// Frequency of light
    Hertz,

    /// Period of light
    Second,
}


/// Each energy quantity should contains three parts: number, prefix and unit.
pub struct Quantity {
    /// Singular float number
    number: f64,

    /// Prefix of the unit.
    /// Example: 'f' in 'fs' (femto second), 'n' in 'ns' (nano second).
    prefix: MetricPrefix,

    /// Energy unit.
    unit: Unit,
}


impl Quantity {
    fn normalize(mut self) -> Self {
        self.normalize_prefix()
            .normalize_unit()
    }

    fn normalize_prefix(mut self) -> Self {
        todo!()
    }

    fn normalize_unit(mut self) -> Self {
        todo!()
    }
}
