use std::sync::OnceLock;
use std::hash::Hash;
use std::collections::HashMap;

use nom::{
    branch::alt,
    bytes::complete::tag,
    character::complete::multispace0 as multispace,
    number::complete::double,
    combinator::map,
    sequence::delimited,
    IResult,
};

use crate::Result;


#[derive(Clone, Copy, PartialEq, PartialOrd, Eq, Ord, Debug, Hash)]
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
    One,

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


macro_rules! _apply {
    ($func:expr, $data:expr) => {
        ( $func($data), )
    };
    ($func:expr $(, $data:expr)+) => {
        ( $( $func($data) ),* )
    };
}

macro_rules! prefix_parser {
    ($prefix:expr $(, $name:expr)+) => (
        map(alt(_apply!(tag, $($name),*)), |_| $prefix)
    );
}


impl MetricPrefix {
    fn parse_prefix(i: &str) -> IResult<&str, MetricPrefix> {
        use MetricPrefix::*;

        let atto  = prefix_parser!(Atto,  "a", "atto",  "Atto");
        let femto = prefix_parser!(Femto, "f", "femto", "Femto");
        let pico  = prefix_parser!(Pico,  "p", "pico",  "Pico");
        let nano  = prefix_parser!(Nano,  "n", "nano",  "Nano");
        let micro = prefix_parser!(Micro, "u", "μ", "mu", "Mu", "micro", "Micro");
        let milli = prefix_parser!(Milli, "m", "milli", "Milli");
        let kilo  = prefix_parser!(Kilo,  "K", "kilo",  "Kilo");
        let mega  = prefix_parser!(Mega,  "M", "mega",  "Mega", "Mi");
        let giga  = prefix_parser!(Giga,  "G", "giga",  "Giga", "Gi");
        let tera  = prefix_parser!(Tera,  "T", "tera",  "Tera", "Ti");
        let peta  = prefix_parser!(Peta,  "P", "peta",  "Peta", "Pi");
        let exa   = prefix_parser!(Exa,   "E", "exa",   "Exa");
        let one   = prefix_parser!(One, "");

        alt((
            atto,
            femto,
            pico,
            nano,
            micro,
            milli,
            kilo,
            mega,
            giga,
            tera,
            peta,
            exa,
            one,
        ))(i)
    }
}


fn get_prefix_scale() -> &'static HashMap<MetricPrefix, f64> {
    static INSTANCE: OnceLock<HashMap<MetricPrefix, f64>> = OnceLock::new();
    &INSTANCE.get_or_init(|| {
        [
            (MetricPrefix::Atto,  1E-18),
            (MetricPrefix::Femto, 1E-15),
            (MetricPrefix::Pico,  1E-12),
            (MetricPrefix::Nano,  1E-9),
            (MetricPrefix::Micro, 1E-6),
            (MetricPrefix::Milli, 1E-3),
            (MetricPrefix::One,   1.0f64),
            (MetricPrefix::Kilo,  1E3),
            (MetricPrefix::Mega,  1E6),
            (MetricPrefix::Giga,  1E9),
            (MetricPrefix::Tera,  1E12),
            (MetricPrefix::Peta,  1E15),
            (MetricPrefix::Exa,   1E18),
        ].iter().cloned().collect()
    })
}


#[derive(Clone, Copy, PartialEq, Eq, PartialOrd, Ord, Debug, Hash)]
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


fn get_ratio_ev_to_other() -> &'static HashMap<Unit, f64> {
    static INSTANCE: OnceLock<HashMap<Unit, f64>> = OnceLock::new();
    &INSTANCE.get_or_init(|| {
        [
            (Unit::ElectronVolt,   1.0f64),
            (Unit::CaloriePerMole, 1.60217733 * 6.0223 * 1E4 / 4.184),
            (Unit::JoulePerMole,   1.60217733 * 6.0223 * 1E4),
            (Unit::Kelvin,         1.160451812E4),
            (Unit::Hartree,        1.0 / 27.2114),
            (Unit::Wavenumber,     8065.73),
            (Unit::Meter,          1.23984193E-6),
            (Unit::Hertz,          2.417989242E14),
            (Unit::Second,         1.0 / 2.417989242E14),
        ].iter().cloned().collect()
    })
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
    pub fn from_str(i: &str) -> Result<Self> {
        todo!();
    }


    pub fn normalize(self) -> Self {
        self.normalize_prefix()
            .normalize_unit()
    }

    fn normalize_prefix(mut self) -> Self {
        let scale = get_prefix_scale()[&self.prefix];
        self.number *= scale;
        self.prefix = MetricPrefix::One;
        self
    }

    // the `prefix` must be `One` before calling this function
    fn normalize_unit(mut self) -> Self {
        use Unit::*;

        //assert_eq!(self.prefix, MetricPrefix::One);
        self = self.normalize_prefix();
        let unit = self.unit;
        let ratio = get_ratio_ev_to_other()[&unit];
        self.number = match unit {
            Meter | Second => ratio / self.number,
            _ => self.number / ratio,
        };
        self.unit = Unit::ElectronVolt;
        self
    }


    fn to_quantity(mut self, unit: Unit) -> Self {
        self.to_normalized_quantity(unit)
            .add_metrix_prefix()
    }

    // the `prefix` must be `One` before calling this function
    fn to_normalized_quantity(mut self, unit: Unit) -> Self {
        use Unit::*;
        //assert_eq!(self.prefix, MetricPrefix::One);
        //assert_eq!(self.unit, Unit::ElectronVolt);
        self = self.normalize();

        self.unit = unit;
        let ratio = get_ratio_ev_to_other()[&unit];
        self.number = match unit {
            Meter | Second => ratio / self.number,
            _ => self.number * ratio,
        };
        self
    }


    fn add_metrix_prefix(mut self) -> Self {
        use MetricPrefix::*;

        //assert_eq!(self.prefix, One);
        self = self.normalize_prefix();
        let number = self.number;
        let prefix = match number {
            x if x <= 1E-18 => Atto,
            _ => Exa,
        };

        self.number /= get_prefix_scale()[&prefix];
        self.prefix  = prefix;
        self
    }
}


fn pnumber(i: &str) -> IResult<&str, f64> {
    delimited(
        multispace,
        map(double, |x| x),
        multispace,
    )(i)
}


fn pprefix(i: &str) -> IResult<&str, MetricPrefix> {
    use MetricPrefix::*;

    alt((
            map(
                alt((
                        tag("a"),
                        tag("atto"),
                        tag("Atto")
                )), |_| Atto
            ),
    ))(i)
}
