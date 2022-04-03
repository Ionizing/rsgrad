use std::{
    io::Write,
    fs,
    path::Path,
    panic::catch_unwind,
};

use serde::{
    Serialize,
    Deserialize,
};
use log::warn;
use anyhow::{
    bail,
    Result,
    Context,
};
use plotly::common::{
    ColorScalePalette,
    color::{
        Color,
        ColorWrapper,
    },
};
use ndarray::Array1;
use structopt::clap::arg_enum;

use crate::types::{
    range_parse,
    index_transform,
};


arg_enum! {
    #[derive(Debug, Clone, Copy, Serialize, Deserialize, PartialEq, Eq)]
    pub enum Axis {
        X,
        Y,
        Z,
    }
}


const NAMED_COLORS: &[&str] = &[
        "aliceblue",            "antiquewhite",     "aqua",             "aquamarine",       "azure",
        "beige",                "bisque",           "black",            "blanchedalmond",   "blue",
        "blueviolet",           "brown",            "burlywood",        "cadetblue",        "chartreuse",
        "chocolate",            "coral",            "cornflowerblue",   "cornsilk",         "crimson",
        "cyan",                 "darkblue",         "darkcyan",         "darkgoldenrod",    "darkgray",
        "darkgrey",             "darkgreen",        "darkkhaki",        "darkmagenta",      "darkolivegreen",
        "darkorange",           "darkorchid",       "darkred",          "darksalmon",       "darkseagreen",
        "darkslateblue",        "darkslategray",    "darkslategrey",    "darkturquoise",    "darkviolet",
        "deeppink",             "deepskyblue",      "dimgray",          "dimgrey",          "dodgerblue",
        "firebrick",            "floralwhite",      "forestgreen",      "fuchsia",          "gainsboro",
        "ghostwhite",           "gold",             "goldenrod",        "gray",             "grey",
        "green",                "greenyellow",      "honeydew",         "hotpink",          "indianred",
        "indigo",               "ivory",            "khaki",            "lavender",         "lavenderblush",
        "lawngreen",            "lemonchiffon",     "lightblue",        "lightcoral",       "lightcyan",
        "lightgoldenrodyellow", "lightgray",        "lightgrey",        "lightgreen",       "lightpink",
        "lightsalmon",          "lightseagreen",    "lightskyblue",     "lightslategray",   "lightslategrey",
        "lightsteelblue",       "lightyellow",      "lime",             "limegreen",        "linen",
        "magenta",              "maroon",           "mediumaquamarine", "mediumblue",       "mediumorchid",
        "mediumpurple",         "mediumseagreen",   "mediumslateblue",  "mediumspringgreen","mediumturquoise",
        "mediumvioletred",      "midnightblue",     "mintcream",        "mistyrose",        "moccasin",
        "navajowhite",          "navy",             "oldlace",          "olive",            "olivedrab",
        "orange",               "orangered",        "orchid",           "palegoldenrod",    "palegreen",
        "paleturquoise",        "palevioletred",    "papayawhip",       "peachpuff",        "peru",
        "pink",                 "plum",             "powderblue",       "purple",           "red",
        "rosybrown",            "royalblue",        "saddlebrown",      "salmon",           "sandybrown",
        "seagreen",             "seashell",         "sienna",           "silver",           "skyblue",
        "slateblue",            "slategray",        "slategrey",        "snow",             "springgreen",
        "steelblue",            "tan",              "teal",             "thistle",          "tomato",
        "turquoise",            "violet",           "wheat",            "white",            "whitesmoke",
        "yellow",               "yellowgreen",      "transparent",
    ];


#[derive(Debug, Clone)]
pub struct CustomColor (ColorWrapper);

impl Color for CustomColor {
    fn to_color(&self) -> ColorWrapper {
        self.0.clone()
    }
}


const PALETTES: &[&str] = &[
    "blackbody",    "bluered",  "blues",    "cividis",  "earth",
    "electric",     "greens",   "greys",    "hot",      "jet",
    "picnic",       "portland", "rainbow",  "rdbu",     "reds",
    "viridis",      "ylgnbu",   "ylorrd",
];


const PALETTES_ENUM: &[ColorScalePalette] = &[
    ColorScalePalette::Blackbody,
    ColorScalePalette::Bluered,
    ColorScalePalette::Blues,
    ColorScalePalette::Cividis,
    ColorScalePalette::Earth,
    ColorScalePalette::Electric,
    ColorScalePalette::Greens,
    ColorScalePalette::Greys,
    ColorScalePalette::Hot,
    ColorScalePalette::Jet,
    ColorScalePalette::Picnic,
    ColorScalePalette::Portland,
    ColorScalePalette::Rainbow,
    ColorScalePalette::RdBu,
    ColorScalePalette::Reds,
    ColorScalePalette::Viridis,
    ColorScalePalette::YlGnBu,
    ColorScalePalette::YlOrRd,
];


pub fn write_array_to_txt(file_name: &(impl AsRef<Path> + ?Sized), ys: Vec<&Array1<f64>>, comment: &str) -> Result<()> {
    let ncol = ys.len();

    let x = ys.get(0).context("At lease two data sets are needed")?;
    let nrow = x.len();

    if nrow == 0 || !ys.iter().all(|y| y.len() == nrow) {
        bail!("[WRT_ARRAY]: input data with zero length or they don't have consistent lengths");
    }

    let mut f = fs::OpenOptions::new()
        .create(true)
        .truncate(true)
        .write(true)
        .open(file_name)?;

    writeln!(f, "# {}", comment.trim())?;

    for irow in 0 .. nrow {
        let mut s = String::with_capacity(8);
        for icol in 0 .. ncol {
            s.push_str(&format!("  {:15.6}", ys[icol][irow]));
        }
        s.push('\n');

        f.write(s.as_bytes())?;
    }
    
    Ok(())
}


#[derive(Clone, Serialize, Deserialize, Debug)]
pub struct RawSelection {
    pub spins:      Option<String>,
    pub kpoints:    Option<String>,
    pub atoms:      Option<String>,
    pub orbits:     Option<String>,
    pub color:      Option<String>,
    pub factor:     Option<f64>,
}


impl RawSelection {
    /// Parse the spin options from configuration, the result is sorted and deduplicated.
    ///
    /// For ispin = 1 system, whatever in, vec![0] out;
    /// For ispin = 2 system, it accepts 'u' 'd' 'up' 'down' 'dn' and uppercased strings;
    /// For NCL system, it accepts 'x' 'y' 'z' 'tot' and uppercased strings.
    pub fn parse_ispins(input: Option<&str>, nspin: usize, is_ncl: bool) -> Result<Vec<usize>> {
        let mut ret = if let Some(spins) = input {
            if spins.trim().is_empty() {
                bail!("[DOS]: No spin component selected.");
            }

            if is_ncl {
                spins.split_whitespace()
                    .map(|x| match x.to_ascii_lowercase().as_ref() {
                        "x"         => Ok(1usize),
                        "y"         => Ok(2usize),
                        "z"         => Ok(3usize),
                        "t" | "tot" => Ok(0usize),
                        _ =>
bail!("[DOS]: Invalid spin component selected: `{}`, available components are `x`, `y`, `z` and `tot`", x)
                    }).collect::<Result<Vec<_>>>()
            } else if nspin == 2 {
                spins.split_whitespace()
                    .map(|x| match x.to_ascii_lowercase().as_ref() {
                        "u" | "up"           => Ok(0usize),
                        "d" | "dn" | "down"  => Ok(1usize),
                        _ =>
bail!("[DOS]: Invalid spin component selected: `{}`, available components are `up` and `down`", x)
                    }).collect::<Result<Vec<_>>>()
            } else {
                warn!("[DOS]: This system is not spin-polarized, only one spin component is available, selected by default.");
                Ok(vec![0usize])
            }
        } else {
            Ok(if is_ncl {
                vec![0usize]  // 'tot' part
            } else if nspin == 2 {
                vec![0usize, 1usize]  // spin up and spin down
            } else {
                vec![0usize]  // only one spin available
            })
        }?;

        ret.sort();
        ret.dedup();
        Ok(ret)
    }

    
    /// Parse the atom index.
    ///
    /// Negative indices are allowed to index from tail. All the indices are sorted and
    /// deduplicated.
    pub fn parse_iatoms(input: Option<&str>, nions: usize) -> Result<Vec<usize>> {
        if let Some(atoms) = input {
            let mut ret = atoms.split_whitespace()
                .map(|x| range_parse(x))
                .collect::<Result<Vec<Vec<i32>>>>()?
                .into_iter()
                .map(|x| index_transform(x, nions).into_iter())
                .flatten()
                .map(|x| (x - 1).rem_euclid(nions))
                .collect::<Vec<usize>>();

            if ret.is_empty() {
                bail!("[DOS]: No atoms selected.");
            }

            ret.sort();
            ret.dedup();
            Ok(ret)
        } else {  // All atoms selected if left blank
            Ok((0 .. nions).collect::<Vec<usize>>())
        }
    }

    /// Parse the kpoint index.
    ///
    /// Negative indices are allowed to index from tail. All the indices are sorted and
    /// deduplicated.
    pub fn parse_ikpoints(input: Option<&str>, nkpoints: usize) -> Result<Vec<usize>> {
        if let Some(kpoints) = input {
            let mut ret = kpoints.split_whitespace()
                .map(|x| range_parse(x))
                .collect::<Result<Vec<Vec<i32>>>>()?
                .into_iter()
                .map(|x| index_transform(x, nkpoints).into_iter())
                .flatten()
                .map(|x| (x - 1).rem_euclid(nkpoints))
                .collect::<Vec<usize>>();

            if ret.is_empty() {
                bail!("[DOS]: No ikpoints selected.");
            }

            ret.sort();
            ret.dedup();
            Ok(ret)
        } else {
            Ok((0 .. nkpoints).collect::<Vec<usize>>())
        }
    }

    /// Parse the orbits' name and convert to orbit index. PROCAR is loaded in advance.
    ///
    /// Indices are sorted and deduplicated.
    pub fn parse_iorbits(input: Option<&str>, nlm: &[String]) -> Result<Vec<usize>> {
        if let Some(orbits) = input {
            if orbits.trim().is_empty() {
                bail!("[DOS]: No orbits selected.");
            }

            let mut ret = orbits.split_whitespace()
                .map(|x| {
                    nlm.iter().position(|x2| x2 == &x)
                        .context(format!("Selected orbit {:?} not available in {:?}", x, &nlm))
                })
            .collect::<Result<Vec<_>>>()?;
            
            ret.sort();
            ret.dedup();

            Ok(ret)
        } else {
            Ok((0 .. nlm.len()).collect::<Vec<usize>>())
        }
    }

    /// Parse the color to this curve.
    pub fn parse_color(input: &str) -> Result<CustomColor> {
        if NAMED_COLORS.contains(&input.to_ascii_lowercase().as_ref()) {
            return Ok(CustomColor(ColorWrapper::S(input.to_owned())))
        } else {
            let ret = catch_unwind(|| {
                input.to_color()
            });

            if ret.is_err() {
                bail!("The input color is neither a named color nor a valid hex code. 
See \"https://developer.mozilla.org/en-US/docs/Web/CSS/color_value for availed named colors.\"");
            } else {
                Ok(CustomColor(ret.unwrap()))
            }
        }
    }

    pub fn parse_colormap(input: &str) -> Result<ColorScalePalette> {
        let input = &input.to_ascii_lowercase();
        if let Some(pos) = PALETTES.iter().position(|x| x == input) {
            Ok(PALETTES_ENUM[pos].to_owned())
        } else {
            bail!("Invlid colormap input, available colormaps: {:?}", PALETTES)
        }
    }
}


#[cfg(test)]
mod test {
    use super::*;

    #[test]
    fn test_parse_ispins() {
        // NCL system
        assert_eq!(RawSelection::parse_ispins(Some("x y z t"),  1, true).unwrap(), vec![0, 1, 2, 3]);
        assert_eq!(RawSelection::parse_ispins(Some("x y z t"),  3, true).unwrap(), vec![0, 1, 2, 3]);
        assert_eq!(RawSelection::parse_ispins(None, 3, true).unwrap(), vec![0]);
        assert!(RawSelection::parse_ispins(Some("x y z s"),     3, true).is_err());
        assert!(RawSelection::parse_ispins(Some(""), 3, true).is_err());
        assert!(RawSelection::parse_ispins(Some("u"),           1, true).is_err());
        assert!(RawSelection::parse_ispins(Some("d"),           1, true).is_err());

        // ISPIN = 2 system
        assert_eq!(RawSelection::parse_ispins(Some("up down u d"),  2, false).unwrap(), vec![0, 1]);
        assert_eq!(RawSelection::parse_ispins(Some("u U up UP"),    2, false).unwrap(), vec![0]);
        assert_eq!(RawSelection::parse_ispins(Some("d D dn DN dN"), 2, false).unwrap(), vec![1]);
        assert_eq!(RawSelection::parse_ispins(None,                 2, false).unwrap(), vec![0, 1]);
        assert!(RawSelection::parse_ispins(Some(""),                2, false).is_err());
        assert!(RawSelection::parse_ispins(Some("x"),               2, false).is_err());
        assert!(RawSelection::parse_ispins(Some("y"),               2, false).is_err());
        assert!(RawSelection::parse_ispins(Some("z"),               2, false).is_err());
        assert!(RawSelection::parse_ispins(Some("t"),               2, false).is_err());

        // ISPIN = 1 system
        assert_eq!(RawSelection::parse_ispins(Some("t"), 1, false).unwrap(), vec![0]);
        assert_eq!(RawSelection::parse_ispins(Some("u"), 1, false).unwrap(), vec![0]);
        assert_eq!(RawSelection::parse_ispins(Some("d"), 1, false).unwrap(), vec![0]);
        assert_eq!(RawSelection::parse_ispins(None,      1, false).unwrap(), vec![0]);
    }


    #[test]
    fn test_parse_iatoms() {
        assert_eq!(RawSelection::parse_iatoms(Some("1..8"), 5).unwrap(), vec![0, 1, 2, 3, 4]);
        assert_eq!(RawSelection::parse_iatoms(Some("-2..-1"), 5).unwrap(), vec![3, 4]);
        assert_eq!(RawSelection::parse_iatoms(None, 5).unwrap(), vec![0, 1, 2, 3, 4]);
        assert!(RawSelection::parse_iatoms(Some("-1..-2"), 5).is_err());
        assert!(RawSelection::parse_iatoms(Some("t"), 5).is_err());
    }

    #[test]
    fn test_parse_ikpoints() {
        assert_eq!(RawSelection::parse_ikpoints(Some("1..8"), 5).unwrap(), vec![0, 1, 2, 3, 4]);
        assert_eq!(RawSelection::parse_ikpoints(Some("-2..-1"), 5).unwrap(), vec![3, 4]);
        assert_eq!(RawSelection::parse_ikpoints(None, 5).unwrap(), vec![0, 1, 2, 3, 4]);
        assert!(RawSelection::parse_ikpoints(Some("-1..-2"), 5).is_err());
        assert!(RawSelection::parse_ikpoints(Some("t"), 5).is_err());
    }

    #[test]
    fn test_parse_iorbits() {
        let nlm = "s     py     pz     px    dxy    dyz    dz2    dxz    dx2" 
            .split_whitespace()
            .map(str::to_string)
            .collect::<Vec<_>>();

        assert_eq!(RawSelection::parse_iorbits(None, &nlm).unwrap(), 
                   (0usize .. 9).collect::<Vec<_>>());
        assert_eq!(RawSelection::parse_iorbits(Some("s dx2"), &nlm).unwrap(), vec![0, 8]);
        assert_eq!(RawSelection::parse_iorbits(Some("s     py     pz     px    dxy    dyz    dz2    dxz    dx2 dx2 s"), &nlm).unwrap(), 
                   (0usize .. 9).collect::<Vec<_>>());
        assert!(RawSelection::parse_iorbits(Some("  \n"), &nlm).is_err());
        assert!(RawSelection::parse_iorbits(Some(" y"), &nlm).is_err());
    }
}
