use std::sync::OnceLock;
use std::hash::Hash;
use std::collections::HashMap;

use nom::{
    branch::alt,
    bytes::complete::tag,
    character::complete::{
        digit1,
        multispace0,
    },
    combinator::{
        opt,
        map,
        recognize,
    },
    //multi::{
        //many1,
    //},
    sequence::{
        delimited,
        tuple,
    },
    IResult,
};

use crate::Result;


#[derive(Clone, Copy, PartialEq, Eq, PartialOrd, Ord, Debug, Hash)]
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

        let atto  = prefix_parser!(Atto,  "atto",  "Atto");
        let femto = prefix_parser!(Femto, "femto", "Femto");
        let pico  = prefix_parser!(Pico,  "pico",  "Pico");
        let nano  = prefix_parser!(Nano,  "nano",  "Nano");
        let micro = prefix_parser!(Micro, "μ",     "micro", "Micro");
        let milli = prefix_parser!(Milli, "milli", "Milli");
        let kilo  = prefix_parser!(Kilo,  "kilo",  "Kilo");
        let mega  = prefix_parser!(Mega,  "mega",  "Mega");
        let giga  = prefix_parser!(Giga,  "giga",  "Giga");
        let tera  = prefix_parser!(Tera,  "tera",  "Tera");
        let peta  = prefix_parser!(Peta,  "peta",  "Peta");
        let exa   = prefix_parser!(Exa,   "exa",   "Exa");

        let atto_abbr  = prefix_parser!(Atto,  "a");
        let femto_abbr = prefix_parser!(Femto, "f");
        let pico_abbr  = prefix_parser!(Pico,  "p");
        let nano_abbr  = prefix_parser!(Nano,  "n");
        let micro_abbr = prefix_parser!(Micro, "Mu", "mu", "u");
        let milli_abbr = prefix_parser!(Milli, "m");
        let kilo_abbr  = prefix_parser!(Kilo,  "K");
        let mega_abbr  = prefix_parser!(Mega,  "Mi", "M");
        let giga_abbr  = prefix_parser!(Giga,  "Gi", "G");
        let tera_abbr  = prefix_parser!(Tera,  "Ti", "T");
        let peta_abbr  = prefix_parser!(Peta,  "Pi", "P");
        let exa_abbr   = prefix_parser!(Exa,   "E");

        //let one   = prefix_parser!(One,   "");


        alt((
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
            )),

            alt((
                atto_abbr,
                femto_abbr,
                pico_abbr,
                nano_abbr,
                micro_abbr,
                milli_abbr,
                kilo_abbr,
                mega_abbr,
                giga_abbr,
                tera_abbr,
                peta_abbr,
                exa_abbr,
            )),

            //one,
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


impl Unit {
    fn parse_unit(i: &str) -> IResult<&str, Unit> {
        use Unit::*;

        let ev         = prefix_parser!(ElectronVolt,   "ElectronVolt");
        let calpmol    = prefix_parser!(CaloriePerMole, "Calorie/mol");
        let jpmol      = prefix_parser!(JoulePerMole,   "Joule/mol");
        let kelvin     = prefix_parser!(Kelvin,         "Kelvin");
        let hartree    = prefix_parser!(Hartree,        "Hartree");
        let wavenumber = prefix_parser!(Wavenumber,     "Cm-1");
        let meter      = prefix_parser!(Meter,          "Meter");
        let hertz      = prefix_parser!(Hertz,          "Hertz");
        let second     = prefix_parser!(Second,         "Second");

        let ev_abbr         = prefix_parser!(ElectronVolt,   "eV");
        let calpmol_abbr    = prefix_parser!(CaloriePerMole, "Cal/mol");
        let jpmol_abbr      = prefix_parser!(JoulePerMole,   "J/mol");
        let kelvin_abbr     = prefix_parser!(Kelvin,         "K");
        let hartree_abbr    = prefix_parser!(Hartree,        "Ha");
        let wavenumber_abbr = prefix_parser!(Wavenumber,     "cm-1");
        let meter_abbr      = prefix_parser!(Meter,          "m");
        let hertz_abbr      = prefix_parser!(Hertz,          "Hz");
        let second_abbr     = prefix_parser!(Second,         "s");

        alt((
            alt((
                ev,
                calpmol,
                jpmol,
                kelvin,
                hartree,
                wavenumber,
                meter,
                hertz,
                second,
            )),
            alt((
                ev_abbr,
                calpmol_abbr,
                jpmol_abbr,
                kelvin_abbr,
                hartree_abbr,
                wavenumber_abbr,
                meter_abbr,
                hertz_abbr,
                second_abbr,
            )),
        ))(i)
    }
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


    fn parse_quantity(i: &str) -> IResult<&str, (f64, MetricPrefix, Unit)> {
        let pprefix = MetricPrefix::parse_prefix;
        let punit   = Unit::parse_unit;

        let with_prefix = tuple((
            double,
            delimited(multispace0, pprefix, multispace0),
            punit,
        ));

        let one = prefix_parser!(MetricPrefix::One, "");
        let without_prefix = tuple((
            double,
            delimited(multispace0, one, multispace0),
            punit,
        ));

        alt((
            with_prefix,
            without_prefix,
        ))(i)
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


fn double(i: &str) -> IResult<&str, f64> {
    fn integral(i: &str) -> IResult<&str, &str> {
        digit1(i)
    }
    fn sign(i: &str) -> IResult<&str, &str> {
        alt((tag("+"), tag("-")))(i)
    }
    fn opt_sign(i: &str) -> IResult<&str, &str> {
        map(opt(sign), |x| x.unwrap_or(""))(i)
    }
    fn fraction(i: &str) -> IResult<&str, &str> {
        recognize(
            tuple((
                tag("."), integral
            ))
        )(i)
    }
    fn exponent(i: &str) -> IResult<&str, &str> {
        recognize(
            tuple((
                alt(( tag("e"), tag("E") )),
                opt_sign,
                integral,
            ))
        )(i)
    }

    map(tuple((
        opt_sign,
        integral,
        map(opt(fraction), |x| x.unwrap_or("") ),
        map(opt(exponent), |x| x.unwrap_or("") ),
    )), |(a, b, c, d)| {
        let s = a.to_string() + b + c + d;
        s.parse::<f64>().unwrap()
    })(i)
}


#[cfg(test)]
mod tests {
    use super::*;


    #[test]
    fn test_parse_prefix() {
        use MetricPrefix::*;

        let parser = MetricPrefix::parse_prefix;
        let cases = vec![
            (Atto,  vec!["atto",  "Atto",  "a"]),
            (Femto, vec!["femto", "Femto", "f"]),
            (Pico,  vec!["pico",  "Pico",  "p"]),
            (Nano,  vec!["nano",  "Nano",  "n"]),
            (Micro, vec!["μ",     "mu",    "Mu", "micro", "Micro", "u"]),
            (Milli, vec!["milli", "Milli", "m"]),
            (Kilo,  vec!["kilo",  "Kilo",  "K"]),
            (Mega,  vec!["mega",  "Mega",  "Mi", "M"]),
            (Giga,  vec!["giga",  "Giga",  "Gi", "G"]),
            (Tera,  vec!["tera",  "Tera",  "Ti", "T"]),
            (Peta,  vec!["peta",  "Peta",  "Pi", "P"]),
            (Exa,   vec!["exa",   "Exa",   "E"]),
        ];

        for (prefix, ss) in cases {
            for s in ss {
                assert_eq!(parser(s), Ok(("", prefix)));
            }
        }
    }

    #[test]
    fn test_parse_unit() {
        use Unit::*;

        let parser = Unit::parse_unit;
        let cases = vec![
            (ElectronVolt,   vec!["ElectronVolt", "eV"]),
            (CaloriePerMole, vec!["Calorie/mol", "Cal/mol"]),
            (JoulePerMole,   vec!["Joule/mol", "J/mol"]),
            (Kelvin,         vec!["Kelvin", "K"]),
            (Hartree,        vec!["Hartree", "Ha"]),
            (Wavenumber,     vec!["Cm-1", "cm-1"]),
            (Meter,          vec!["Meter", "m"]),
            (Hertz,          vec!["Hertz", "Hz"]),
            (Second,         vec!["Second", "s"]),
        ];

        for (unit, ss) in cases {
            for s in ss {
                assert_eq!(parser(s), Ok(("", unit)));
            }
        }
    }

    #[test]
    fn test_parse_quantity() {
        use MetricPrefix::*;
        use Unit::*;

        let parser = Quantity::parse_quantity;

        let prefix_cases = vec![
            (Atto,  vec!["atto",  "Atto",  "a"]),
            (Femto, vec!["femto", "Femto", "f"]),
            (Pico,  vec!["pico",  "Pico",  "p"]),
            (Nano,  vec!["nano",  "Nano",  "n"]),
            (Micro, vec!["μ",     "mu",    "Mu", "micro", "Micro", "u"]),
            (Milli, vec!["milli", "Milli", "m"]),
            (Kilo,  vec!["kilo",  "Kilo",  "K"]),
            (Mega,  vec!["mega",  "Mega",  "Mi", "M"]),
            (Giga,  vec!["giga",  "Giga",  "Gi", "G"]),
            (Tera,  vec!["tera",  "Tera",  "Ti", "T"]),
            (Peta,  vec!["peta",  "Peta",  "Pi", "P"]),
            (Exa,   vec!["exa",   "Exa",   "E"]),
        ];

        let unit_cases = vec![
            (ElectronVolt,   vec!["ElectronVolt", "eV"]),
            (CaloriePerMole, vec!["Calorie/mol", "Cal/mol"]),
            (JoulePerMole,   vec!["Joule/mol", "J/mol"]),
            (Kelvin,         vec!["Kelvin", "K"]),
            (Hartree,        vec!["Hartree", "Ha"]),
            (Wavenumber,     vec!["Cm-1", "cm-1"]),
            (Meter,          vec!["Meter", "m"]),
            (Hertz,          vec!["Hertz", "Hz"]),
            (Second,         vec!["Second", "s"]),
        ];

        for (prefix, ssprefix) in &prefix_cases {
            for (unit, ssunit) in &unit_cases {
                for sunit in ssunit {
                    for sprefix in ssprefix {
                        let s = vec!["1.0", sprefix, sunit].join(" ");
                        assert_eq!(parser(&s), Ok(("", (1.0f64, *prefix, *unit))), "{}", s);

                        let s = vec!["0.05E2", sprefix, sunit].join(" ");
                        assert_eq!(parser(&s), Ok(("", (5.0f64, *prefix, *unit))), "{}", s);

                        let s = vec!["1.0", sprefix, sunit].join("");
                        assert_eq!(parser(&s), Ok(("", (1.0f64, *prefix, *unit))), "{}", s);

                        let s = vec!["0.05E2", sprefix, sunit].join("");
                        assert_eq!(parser(&s), Ok(("", (5.0f64, *prefix, *unit))), "{}", s);
                    }

                    let s = vec!["1.0", sunit].join(" ");
                    assert_eq!(parser(&s), Ok(("", (1.0f64, One, *unit))), "{}", s);

                    let s = vec!["0.05E2", sunit].join(" ");
                    assert_eq!(parser(&s), Ok(("", (5.0f64, One, *unit))), "{}", s);

                    let s = vec!["1.0", sunit].join("");
                    assert_eq!(parser(&s), Ok(("", (1.0f64, One, *unit))), "{}", s);

                    let s = vec!["0.05E2", sunit].join("");
                    assert_eq!(parser(&s), Ok(("", (5.0f64, One, *unit))), "{}", s);
                }
            }
        }

    }
}
