use std::{
    io::{
        self,
        BufRead,
        Seek,
        SeekFrom,
        BufReader,
    },
    path::Path,
    fs::File,
};

use regex::Regex;
use ndarray::{
    Array1,
    Array3,
};
use anyhow::{
    Context,
    bail,
};
use rayon::prelude::*;

use crate::{
    traits::Result,
    Poscar,
};


/// Main struct of volumetric data
///
/// # CHGCAR
///
/// This file contains the lattice vectors, atomic coordinates, the total charge density multiplied
/// by the volume `rho(r) * V_cell` on the fine FFT-grid `(NG(X,Y,Z)F)`, and the PAW
/// one-center occupancies. CHGCAR can be used to restart VASP from an existing charge density.
///
/// ## Structure of CHGCAR
///
/// Here is a 'pseudo' CHGCAR file content
///
/// ```text
/// unknown system                          \
/// 1.00000000000000                        |
/// 2.969072   -0.000523   -0.000907        |
/// -0.987305    2.800110    0.000907       |
/// -0.987305   -1.402326    2.423654       |-> positions of atom in POSCAR format
/// Li                                      |
/// 1                                       |
/// Direct                                  |
/// 0.000000  0.000000  0.000000            /
///
/// 2    3    4                             |-> number of grids in x, y, z directions.
/// 0.44 0.44 0.46 0.48 0.52   \
/// 0.56 0.60 0.66 0.73 0.80   |
/// 0.88 0.94 0.10 0.10 0.10   |-> Total charge density
/// 0.10 0.10 0.10 0.10 0.10   |
/// 0.10 0.10 0.10 0.10        /
/// augmentation occupancies 1 15  \
/// 0.27 -0.33  0.00  0.00  0.00   |
/// 0.10  0.00  0.00  0.00  0.39   |
/// 0.58 -0.72 -0.36  0.10 -0.20   |-> PAW augmentation data
/// augmentation occupancies 2 15  |
/// 0.27 -0.33  0.00  0.00  0.00   |
/// 0.10  0.00  0.00  0.00  0.39   |
/// 0.58 -0.72 -0.36  0.10 -0.20   /
/// 2    3    4                             |-> number of grids in x, y, z directions.
/// 0.44 0.44 0.46 0.48 0.52   \
/// 0.56 0.60 0.66 0.73 0.80   |    rho(up) - rho(dn) in ISPIN=2 system
/// 0.88 0.94 0.10 0.10 0.10   | -> rho_x in non collinear system
/// 0.10 0.10 0.10 0.10 0.10   |    NO THIS PART IN ISPIN=1 SYSTEM
/// 0.10 0.10 0.10 0.12        /
/// augmentation occupancies 1 15  \
/// 0.27 -0.33  0.00  0.00  0.00   |
/// 0.10  0.00  0.00  0.00  0.39   |
/// 0.58 -0.72 -0.36  0.10 -0.20   |
/// augmentation occupancies 2 15  |-> PAW augmentation data
/// 0.27 -0.33  0.00  0.00  0.00   |
/// 0.10  0.00  0.00  0.00  0.39   |
/// 0.58 -0.72 -0.36  0.10 -0.00   /
/// <-- If this is an SOC system, another TWO charge density difference parts should be in the following -->
/// <-- NGX NGY NGZ -->  rho_y
/// <-- GRID DATA -->
/// <-- NGX NGY NGZ -->  rho_z
/// <-- GRID DATA -->
/// ```
///
/// ## Structure of PARCHG/CHG
///
/// Similar to the structure of CHGCAR, but without augmentation parts.
///
/// PARCHG is the partial charge density which only takes the charge density of
/// energy/band/kpoint specified electron states.
///
/// Also, CHG stores the total charge density of all the electrons below fermi level in all kpoint,
/// all bands.
///
struct ChargeDensity {
    pub pos:        Poscar,
    pub ngrid:      [usize; 3],
    pub chg:        Vec<Array3<f64>>,
    pub aug:        Vec<String>,
}


impl ChargeDensity {
    pub fn from_str(txt: &str) -> Result<Self> {
        let separate_pos = Regex::new(r"(?m)^\s*$").unwrap()
            .find(txt)
            .context("[CHG]: This file has no empty line to separate position data and grid data.")?
            .start();

        let pos = Self::read_poscar(txt)?;

        let chg_starts = Regex::new(r"(?m)^\s+\d+\s+\d+\s+\d+\s*$").unwrap()
            .find_iter(&txt[separate_pos ..])
            .map(|m| m.start() + separate_pos)
            .collect::<Vec<usize>>();

        if chg_starts.len() == 0 {
            bail!("[CHG]: This file has no grid size data.");
        }

        let chg = chg_starts.par_iter()
            .map(|p| Self::read_chg(&txt[*p ..]))
            .collect::<Result<Vec<Array3<f64>>>>()?;

        let ngrid = {
            let s = chg[0].shape();
            [s[0], s[1], s[2]]
        };

        let aug = chg_starts.par_iter()
            .map(|p| Self::read_raw_aug(&txt[*p ..]))
            .collect::<Option<Vec<String>>>().unwrap_or(vec![]);

        Ok(Self {
            pos,
            ngrid,
            chg,
            aug,
        })
    }


    /// Read CHGCAR header to get POSCAR info
    fn read_poscar(txt: &str) -> Result<Poscar> {
        Poscar::from_str(txt)
    }

    /// This function reads the grid data.  input str should be like:
    ///
    ///  2 3 4
    ///  ... ...
    ///  ... ...
    ///  augmentation ....
    ///  ...
    ///  ...
    ///
    ///  The augmentation part is dropped
    fn read_chg(txt: &str) -> Result<Array3<f64>> {
        let mat = Regex::new(r"(?m)^\s+\d+\s+\d+\s+\d+\s*$").unwrap().find(txt).unwrap();
        let ngrid = {
            let mut v = vec![0usize; 3];
            for (i, s) in txt[mat.start() .. mat.end()].split_whitespace()
                    .take(3).enumerate() {
                v[i] = s.parse::<usize>()
                    .context(format!("[CHG]: Cannot parse {} into float number", s))?;
            }
            (v[0], v[1], v[2])
        };
        
        let start_pos = mat.end();
        let aug_pos = txt.find("augmentation").unwrap_or(txt.len());

        let chg = txt[start_pos .. aug_pos]
            .split_whitespace()
            .map(|s| s.parse::<f64>().expect(&format!("[CHG]: Cannot parse {} into float number", s)))
            .collect::<Array1<f64>>()
            .into_shape(ngrid)?;

        Ok(chg)
    }

    fn read_raw_aug(txt: &str) -> Option<String> {
        let start_pos = txt.find("augmentation")?;
        let end_pos = Regex::new(r"(?m)^\s+\d+\s+\d+\s+\d+\s*$")  // find NGX NGY NGZ
            .unwrap()
            .find(&txt[start_pos ..])
            .map(|x| x.start() + start_pos)
            .unwrap_or(txt.len());

        Some( txt[start_pos .. end_pos].to_string() )
    }
}



#[cfg(test)]
mod test {
    use super::*;

    const SAMPLE: &'static str = "\
unknown system
   1.00000000000000
     2.969072   -0.000523   -0.000907
    -0.987305    2.800110    0.000907
    -0.987305   -1.402326    2.423654
   Li
     1
Direct
  0.000000  0.000000  0.000000

    2    3    4
 0.44062142953E+00 0.44635237036E+00 0.46294638829E+00 0.48881056285E+00 0.52211506729E+00
 0.56203432815E+00 0.60956087775E+00 0.66672131696E+00 0.73417916031E+00 0.80884817972E+00
 0.88351172791E+00 0.94912993844E+00 0.10000382501E+01 0.10353398391E+01 0.10568153616E+01
 0.10677009023E+01 0.10709392990E+01 0.10677009023E+01 0.10568153616E+01 0.10353398391E+01
 0.10677009023E+01 0.10709392990E+01 0.10677009023E+01 0.10568153616E+01
augmentation occupancies 1 15
  0.2743786E+00 -0.3307158E-01  0.0000000E+00  0.0000000E+00  0.0000000E+00
  0.1033253E-02  0.0000000E+00  0.0000000E+00  0.0000000E+00  0.3964234E-01
  0.5875445E-05 -0.7209739E-05 -0.3625569E-05  0.1019266E-04 -0.2068344E-05
augmentation occupancies 2 15
  0.2743786E+00 -0.3307158E-01  0.0000000E+00  0.0000000E+00  0.0000000E+00
  0.1033253E-02  0.0000000E+00  0.0000000E+00  0.0000000E+00  0.3964234E-01
  0.5875445E-05 -0.7209739E-05 -0.3625569E-05  0.1019266E-04 -0.2068344E-05
    2    3    4
 0.44062142953E+00 0.44635237036E+00 0.46294638829E+00 0.48881056285E+00 0.52211506729E+00
 0.56203432815E+00 0.60956087775E+00 0.66672131696E+00 0.73417916031E+00 0.80884817972E+00
 0.88351172791E+00 0.94912993844E+00 0.10000382501E+01 0.10353398391E+01 0.10568153616E+01
 0.10677009023E+01 0.10709392990E+01 0.10677009023E+01 0.10568153616E+01 0.10353398391E+01
 0.10677009023E+01 0.10709392990E+01 0.10677009023E+01 0.12668153616E+01
augmentation occupancies 1 15
  0.2743786E+00 -0.3307158E-01  0.0000000E+00  0.0000000E+00  0.0000000E+00
  0.1033253E-02  0.0000000E+00  0.0000000E+00  0.0000000E+00  0.3964234E-01
  0.5875445E-05 -0.7209739E-05 -0.3625569E-05  0.1019266E-04 -0.2038144E-05
augmentation occupancies 2 15
  0.2743786E+00 -0.3307158E-01  0.0000000E+00  0.0000000E+00  0.0000000E+00
  0.1033253E-02  0.0000000E+00  0.0000000E+00  0.0000000E+00  0.3964234E-01
  0.5875445E-05 -0.7209739E-05 -0.3625569E-05  0.1019266E-04 -0.0038244E-05
";
    
    #[test]
    fn test_read_poscar() {
        let pos = ChargeDensity::read_poscar(SAMPLE).unwrap();
    }

    #[test]
    fn test_read_chg() {
        let chg = ChargeDensity::read_chg(&SAMPLE[200..]).unwrap();
        assert_eq!(chg.shape(), &[2usize, 3, 4]);
        assert_eq!(chg[[1, 2, 3]], 0.10568153616E+01);
    }

    #[test]
    fn test_read_raw_aug() {
        let aug = ChargeDensity::read_raw_aug(SAMPLE);
        assert!(aug.is_some());
    }

    #[test]
    fn test_from_str() {
        let chg = ChargeDensity::from_str(SAMPLE).unwrap();
        assert_eq!(chg.ngrid, [2, 3, 4]);
        assert_eq!(chg.chg.len(), 2);
        assert_eq!(chg.aug.len(), 2);
    }
}
