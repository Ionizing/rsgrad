use std::{
    path::Path,
    fs,
    fmt,
    ops::{
        Add,
        Sub,
    },
};

use regex::Regex;
use ndarray::{
    Array3,
    ShapeBuilder,
};
use anyhow::{
    Context,
    bail,
};
use rayon::prelude::*;

use crate::{
    Result,
    Poscar,
};


#[derive(Clone, Copy, Debug, PartialEq, Eq)]
pub enum ChargeType {
    Chgcar,
    Locpot,
}



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
#[derive(Clone, Debug)]
pub struct ChargeDensity {
    pub chgtype:    ChargeType,
    pub pos:        Poscar,
    pub ngrid:      [usize; 3],
    pub chg:        Vec<Array3<f64>>,
    pub aug:        Vec<String>,
}


impl ChargeDensity {
    /// Read CHGCAR like volumetric data from file.
    pub fn from_file(path: &(impl AsRef<Path> + ?Sized), chgtype: ChargeType) -> Result<Self> {
        let txt = fs::read_to_string(path)?;
        Self::from_str(&txt, chgtype)
    }


    /// Parse CHGCAR like volumetric data from string.
    pub fn from_str(txt: &str, chgtype: ChargeType) -> Result<Self> {
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

        let mut chg = chg_starts.par_iter()
            .map(|p| Self::read_chg(&txt[*p ..]))
            .collect::<Result<Vec<Array3<f64>>>>()?;

        match chgtype {
            ChargeType::Chgcar => {
                let vol = pos.get_volume();
                chg.par_iter_mut()
                    .for_each(|charge| *charge /= vol);
            },
            ChargeType::Locpot => { }
        }

        let ngrid = {
            let s = chg[0].shape();
            [s[0], s[1], s[2]]
        };

        let aug = chg_starts.par_iter()
            .map(|p| Self::read_raw_aug(&txt[*p ..]))
            .collect::<Option<Vec<String>>>().unwrap_or(vec![]);

        if !aug.is_empty() && chg.len() != aug.len() {
            bail!("[CHG]: Augmentation data sets' count not consistent with chage densities' : {} != {}", aug.len(), chg.len());
        }

        Ok(Self {
            chgtype,
            pos,
            ngrid,
            chg,
            aug,
        })
    }


    pub fn to_file(&self, path: &(impl AsRef<Path> + ?Sized)) -> Result<()> {
        fs::write(path, self.to_string())?;
        Ok(())
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
        let _regex = Regex::new(r"(?m)^\s+\d+\s+\d+\s+\d+\s*$").unwrap();
        let mat = _regex.find(txt).unwrap();
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
        let end_pos = txt.find("augmentation")  // end_pos should be the start of augmentation or next NGXF NGYF NGZF or end of txt
            .or_else(|| {
                let mat = _regex.find(&txt[start_pos ..])?;
                Some(mat.start() + start_pos)
            })
            .unwrap_or(txt.len());

        let chg_vec = txt[start_pos .. end_pos]
            .split_whitespace()
            .map(|s| s.parse::<f64>().expect(&format!("[CHG]: Cannot parse {} into float number", s)))
            .collect::<Vec<f64>>();
        let chg = Array3::from_shape_vec(ngrid.f(), chg_vec)?;

        Ok(chg)
    }


    /// All augmentation data is extracted without parsing. The tail linebreaks are preserved.
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


impl fmt::Display for ChargeDensity {
    fn fmt(&self, f: &mut fmt::Formatter<'_>) -> fmt::Result {

        // Flatten chg struct
        let chg = match self.chgtype {
            ChargeType::Chgcar => {
                let chg = self.chg.clone();
                let vol = self.pos.get_volume();
                chg.into_par_iter()
                    .map(|mut charge| {
                        charge *= vol;
                        charge.into_raw_vec()
                    })
                    .collect::<Vec<_>>()
            },
            ChargeType::Locpot => {
                self.chg.par_iter()
                    .map(|charge| charge.clone().into_raw_vec())
                    .collect::<Vec<_>>()
            }
        };

        writeln!(f, "{}", self.pos.to_formatter())?;  // contains the empty line for seperation already

        let empty_aug = vec!["".to_string()];
        let aug = if self.aug.len() != 0 {
            &self.aug
        } else {
            &empty_aug
        };

        for (c, a) in chg.iter().zip(aug.iter().cycle()) {
            // write NGFX NGFY NGFZ
            writeln!(f, " {:5} {:5} {:5}", self.ngrid[0], self.ngrid[1], self.ngrid[2])?;
            
            // write the grid data
            for chunk in c.chunks(5) {
                for ch in chunk.iter() {
                    write!(f, " {:17.11E}", ch)?;
                }
                writeln!(f)?;
            }
            
            f.write_str(a)?;
        }

        Ok(())
    }
}


impl Add for ChargeDensity {
    type Output=Result<Self>;

    /// Add two provided ChargeDensity object, they should have same cell and grid size.
    /// The augmentation part will be dropped.
    fn add(mut self, mut other:Self) -> Self::Output {
        self.pos = self.pos.normalize();
        other.pos = other.pos.normalize();

        if self.chgtype != other.chgtype {
            bail!("[CHG_ADD]: Different type of charge densities are provided: {:?} != {:?}",
                  self.chgtype, other.chgtype);
        }

        if self.pos.cell != other.pos.cell {
            bail!("[CHG_ADD]: Cannot add two system's charge densities within different lattices: \n{:#?}\n != \n{:#?}\n",
                  self.pos.cell, other.pos.cell);
        }

        if self.ngrid != other.ngrid {
            bail!("[CHG_ADD]: Cannot add two system's charge densities within different grids: \n{:#?}\n != \n{:#?}\n",
                  self.ngrid, other.ngrid);
        }

        if self.chg.len() == 0 || self.chg.len() != other.chg.len() {
            bail!("[CHG_ADD]: No charge densitie found or two systems' charge densitie set counts not match: {} != {}",
                  self.chg.len(), other.chg.len());
        }

        if self.pos.constraints.is_some() != other.pos.constraints.is_some() {
            bail!("[CHG_ADD]: Not all provided charge densities have constraints");
        }

        // Construct POSCAR
        let pos = {
            let comment         = "Added charge density. Produced by rsgrad".to_string();
            let scale           = 1.0f64;
            let cell            = self.pos.cell;
            let ion_types       = self.pos.ion_types.into_iter()
                .chain(other.pos.ion_types.into_iter())
                .collect::<Vec<_>>();
            let ions_per_type   = self.pos.ions_per_type.into_iter()
                .chain(other.pos.ions_per_type.into_iter())
                .collect::<Vec<_>>();
            
            let pos_cart        = self.pos.pos_cart.into_iter()
                .chain(other.pos.pos_cart.into_iter())
                .collect::<Vec<_>>();
            let pos_frac        = self.pos.pos_frac.into_iter()
                .chain(other.pos.pos_frac.into_iter())
                .collect::<Vec<_>>();
            let constraints     = self.pos.constraints
                .zip(other.pos.constraints)
                .map(|(x, y)| {
                    x.into_iter()
                     .chain(y.into_iter()).collect::<Vec<_>>()
                });

            if ions_per_type.iter().sum::<i32>() != pos_cart.len() as i32 || 
                pos_cart.len() != pos_frac.len() || 
                (constraints.is_some() && constraints.as_ref().unwrap().len() != pos_cart.len()) {
                bail!("[CHG_ADD]: Invalid POSCAR generated. This is a bug, please open an issue on github and the author will help you solve it.");
            } 

            Poscar {
                comment,
                scale,
                cell,
                ion_types,
                ions_per_type,
                pos_cart,
                pos_frac,
                constraints,
            }
        };

        let chgtype = self.chgtype;
        let ngrid = self.ngrid;
        let chg = self.chg.into_iter().zip(other.chg.into_iter())
            .map(|(x, y)| x + y)
            .collect::<Vec<_>>();

        Ok(Self {
            chgtype,
            pos,
            ngrid,
            chg,
            aug: vec![],
        })
    }
}


impl Sub for ChargeDensity {
    type Output=Result<Self>;

    /// Subtract self's charge density from other's.
    fn sub(mut self, mut other: Self) -> Self::Output {
        self.pos = self.pos.normalize();
        other.pos = other.pos.normalize();

        if self.chgtype != other.chgtype {
            bail!("[CHG_SUB]: Different type of charge densities are provided: {:?} != {:?}",
                  self.chgtype, other.chgtype);
        }

        if self.pos.cell != other.pos.cell {
            bail!("[CHG_SUB]: Cannot subtract two system's charge densities within different lattices: \n{:#?}\n != \n{:#?}\n",
                  self.pos.cell, other.pos.cell);
        }

        if self.ngrid != other.ngrid {
            bail!("[CHG_SUB]: Cannot subtract two system's charge densities within different grids: \n{:#?}\n != \n{:#?}\n",
                  self.ngrid, other.ngrid);
        }

        if self.chg.len() == 0 || self.chg.len() != other.chg.len() {
            bail!("[CHG_SUB]: No charge densitie found or two systems' charge densitie set counts not match: {} != {}",
                  self.chg.len(), other.chg.len());
        }

        if self.pos.constraints.is_some() != other.pos.constraints.is_some() {
            bail!("[CHG_SUB]: Not all provided charge densities have constraints");
        }

        // Construct POSCAR, here we choose the system with more atoms as the final structure.
        let pos = if self.pos.get_natoms() > other.pos.get_natoms() {
            self.pos
        } else {
            other.pos
        };

        let chgtype = self.chgtype;
        let ngrid = self.ngrid;
        let chg = self.chg.into_iter().zip(other.chg.into_iter())
            .map(|(x, y)| x - y)
            .collect::<Vec<_>>();

        Ok(Self{
            chgtype,
            pos,
            ngrid,
            chg,
            aug: vec![],
        })
    }
}


#[cfg(test)]
mod test {
    use super::*;

    const SAMPLE_CHGCAR: &'static str = "\
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
    

    const SAMPLE_CHG: &'static str = "\
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
    2    3    4
 0.44062142953E+00 0.44635237036E+00 0.46294638829E+00 0.48881056285E+00 0.52211506729E+00
 0.56203432815E+00 0.60956087775E+00 0.66672131696E+00 0.73417916031E+00 0.80884817972E+00
 0.88351172791E+00 0.94912993844E+00 0.10000382501E+01 0.10353398391E+01 0.10568153616E+01
 0.10677009023E+01 0.10709392990E+01 0.10677009023E+01 0.10568153616E+01 0.10353398391E+01
 0.10677009023E+01 0.10709392990E+01 0.10677009023E+01 0.12668153616E+01
";


    #[test]
    fn test_read_poscar() {
        ChargeDensity::read_poscar(SAMPLE_CHGCAR).unwrap();
        ChargeDensity::read_poscar(SAMPLE_CHG).unwrap();
    }

    #[test]
    fn test_read_chg() {
        let chg = ChargeDensity::read_chg(&SAMPLE_CHGCAR[200..]).unwrap();
        assert_eq!(chg.shape(), &[2usize, 3, 4]);
        assert_eq!(chg[[0, 1, 0]], 0.46294638829E+00);
        assert_eq!(chg[[1, 2, 3]], 0.10568153616E+01);
        assert_eq!(chg[[0, 0, 0]], 0.44062142953E+00);

        ChargeDensity::read_chg(&SAMPLE_CHG[200..]).unwrap();
    }

    #[test]
    fn test_read_raw_aug() {
        let aug = ChargeDensity::read_raw_aug(SAMPLE_CHGCAR);
        assert!(aug.is_some());
        let aug = aug.unwrap();
        assert!(aug.starts_with("augmentation occupancies 1 15"));
        assert!(aug.ends_with("-0.2068344E-05\n"));

        let aug = ChargeDensity::read_raw_aug(SAMPLE_CHG);
        assert!(aug.is_none());
    }

    #[test]
    fn test_from_str() {
        let chg = ChargeDensity::from_str(SAMPLE_CHGCAR, ChargeType::Chgcar).unwrap();
        assert_eq!(chg.ngrid, [2, 3, 4]);
        assert_eq!(chg.chg.len(), 2);
        assert_eq!(chg.aug.len(), 2);

        let chg = ChargeDensity::from_str(SAMPLE_CHG, ChargeType::Chgcar).unwrap();
        assert_eq!(chg.ngrid, [2, 3, 4]);
        assert_eq!(chg.chg.len(), 2);
        assert_eq!(chg.aug.len(), 0);
    }

    #[test]
    #[ignore]
    fn test_from_file() {
        ChargeDensity::from_file("/Users/ionizing/tmp/CHGCAR", ChargeType::Chgcar).unwrap();
    }

    #[test]
    fn test_format() {
        let format_expect = "\
unknown system
 1.0000000
       2.969072000      -0.000523000      -0.000907000
      -0.987305000       2.800110000       0.000907000
      -0.987305000      -1.402326000       2.423654000
     Li
      1
Direct
      0.0000000000      0.0000000000      0.0000000000 !     Li-001    1

     2     3     4
  4.40621429530E-1  4.46352370360E-1  4.62946388290E-1  4.88810562850E-1  5.22115067290E-1
  5.62034328150E-1  6.09560877750E-1  6.66721316960E-1  7.34179160310E-1  8.08848179720E-1
  8.83511727910E-1  9.49129938440E-1   1.00003825010E0   1.03533983910E0   1.05681536160E0
   1.06770090230E0   1.07093929900E0   1.06770090230E0   1.05681536160E0   1.03533983910E0
   1.06770090230E0   1.07093929900E0   1.06770090230E0   1.05681536160E0
augmentation occupancies 1 15
  0.2743786E+00 -0.3307158E-01  0.0000000E+00  0.0000000E+00  0.0000000E+00
  0.1033253E-02  0.0000000E+00  0.0000000E+00  0.0000000E+00  0.3964234E-01
  0.5875445E-05 -0.7209739E-05 -0.3625569E-05  0.1019266E-04 -0.2068344E-05
augmentation occupancies 2 15
  0.2743786E+00 -0.3307158E-01  0.0000000E+00  0.0000000E+00  0.0000000E+00
  0.1033253E-02  0.0000000E+00  0.0000000E+00  0.0000000E+00  0.3964234E-01
  0.5875445E-05 -0.7209739E-05 -0.3625569E-05  0.1019266E-04 -0.2068344E-05
     2     3     4
  4.40621429530E-1  4.46352370360E-1  4.62946388290E-1  4.88810562850E-1  5.22115067290E-1
  5.62034328150E-1  6.09560877750E-1  6.66721316960E-1  7.34179160310E-1  8.08848179720E-1
  8.83511727910E-1  9.49129938440E-1   1.00003825010E0   1.03533983910E0   1.05681536160E0
   1.06770090230E0   1.07093929900E0   1.06770090230E0   1.05681536160E0   1.03533983910E0
   1.06770090230E0   1.07093929900E0   1.06770090230E0   1.26681536160E0
augmentation occupancies 1 15
  0.2743786E+00 -0.3307158E-01  0.0000000E+00  0.0000000E+00  0.0000000E+00
  0.1033253E-02  0.0000000E+00  0.0000000E+00  0.0000000E+00  0.3964234E-01
  0.5875445E-05 -0.7209739E-05 -0.3625569E-05  0.1019266E-04 -0.2038144E-05
augmentation occupancies 2 15
  0.2743786E+00 -0.3307158E-01  0.0000000E+00  0.0000000E+00  0.0000000E+00
  0.1033253E-02  0.0000000E+00  0.0000000E+00  0.0000000E+00  0.3964234E-01
  0.5875445E-05 -0.7209739E-05 -0.3625569E-05  0.1019266E-04 -0.0038244E-05
";

        let chg = ChargeDensity::from_str(SAMPLE_CHGCAR, ChargeType::Chgcar).unwrap();
        assert_eq!(chg.to_string(), format_expect);
        ChargeDensity::from_str(format_expect, ChargeType::Chgcar).unwrap();
        
        let format_expect = "\
unknown system
 1.0000000
       2.969072000      -0.000523000      -0.000907000
      -0.987305000       2.800110000       0.000907000
      -0.987305000      -1.402326000       2.423654000
     Li
      1
Direct
      0.0000000000      0.0000000000      0.0000000000 !     Li-001    1

     2     3     4
  4.40621429530E-1  4.46352370360E-1  4.62946388290E-1  4.88810562850E-1  5.22115067290E-1
  5.62034328150E-1  6.09560877750E-1  6.66721316960E-1  7.34179160310E-1  8.08848179720E-1
  8.83511727910E-1  9.49129938440E-1   1.00003825010E0   1.03533983910E0   1.05681536160E0
   1.06770090230E0   1.07093929900E0   1.06770090230E0   1.05681536160E0   1.03533983910E0
   1.06770090230E0   1.07093929900E0   1.06770090230E0   1.05681536160E0
     2     3     4
  4.40621429530E-1  4.46352370360E-1  4.62946388290E-1  4.88810562850E-1  5.22115067290E-1
  5.62034328150E-1  6.09560877750E-1  6.66721316960E-1  7.34179160310E-1  8.08848179720E-1
  8.83511727910E-1  9.49129938440E-1   1.00003825010E0   1.03533983910E0   1.05681536160E0
   1.06770090230E0   1.07093929900E0   1.06770090230E0   1.05681536160E0   1.03533983910E0
   1.06770090230E0   1.07093929900E0   1.06770090230E0   1.26681536160E0
";

        let chg = ChargeDensity::from_str(SAMPLE_CHG, ChargeType::Chgcar).unwrap();
        assert_eq!(chg.to_string(), format_expect);
        ChargeDensity::from_str(format_expect, ChargeType::Chgcar).unwrap();
    }

    #[test]
    fn test_chg_add() {
        let chg1 = ChargeDensity::from_str(SAMPLE_CHGCAR, ChargeType::Chgcar).unwrap();
        let chg2 = ChargeDensity::from_str(SAMPLE_CHG, ChargeType::Chgcar).unwrap();

        let chgdat = chg1.chg[0].clone() * 2.0;

        let chg3 = (chg1 + chg2).unwrap();
        assert_eq!(chg3.pos.get_natoms(), 2);
        assert_eq!(chg3.chg[0], chgdat);
    }

    #[test]
    #[should_panic]
    fn test_chg_add_failed() {
        let chg1 = ChargeDensity::from_str(SAMPLE_CHGCAR, ChargeType::Chgcar).unwrap();
        let chg2 = ChargeDensity::from_str(SAMPLE_CHG, ChargeType::Locpot).unwrap();
        (chg1 + chg2).unwrap();
    }

    #[test]
    fn test_chg_sub() {
        let chg1 = ChargeDensity::from_str(SAMPLE_CHGCAR, ChargeType::Chgcar).unwrap();
        let chg2 = ChargeDensity::from_str(SAMPLE_CHG, ChargeType::Chgcar).unwrap();
        
        let chg3 = (chg1 - chg2).unwrap();
        assert!(chg3.chg[0].iter().all(|x| *x == 0.0f64));
    }

    #[test]
    #[should_panic]
    fn test_chg_sub_failed() {
        let chg1 = ChargeDensity::from_str(SAMPLE_CHGCAR, ChargeType::Chgcar).unwrap();
        let chg2 = ChargeDensity::from_str(SAMPLE_CHG, ChargeType::Locpot).unwrap();
        (chg1 - chg2).unwrap();
    }
}
