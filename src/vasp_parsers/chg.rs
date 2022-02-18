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

use ndarray::Array3;
use anyhow::bail;

use crate::{
    traits::Result,
    Poscar,
};

struct ChgBase {
    pos:        Poscar,
    ngrid:      [usize; 3],
    chg:        Array3<f64>,

    // Optional parts
    aug:        Option<String>,
    chgdiff:    Vec<Array3<f64>>,
    augdiff:    Vec<String>,
}


impl ChgBase {
    fn _read_poscar(file: &mut (impl BufRead+Seek)) -> Result<Poscar> {
        let mut buf = String::new();
        while let Ok(n) = file.read_line(&mut buf) {
            if n + 1 == buf.len() - buf.trim_end().len() {
                break;
            }
        }

        Poscar::from_str(&buf)
    }

    fn _read_chg(file: &mut (impl BufRead+Seek)) -> Result<Array3<f64>> {
        todo!();
    }

    fn _read_raw_aug(file: &mut (impl BufRead+Seek)) -> Result<String> {
        todo!();
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
        let mut stream = io::Cursor::new(SAMPLE.as_bytes());
        ChgBase::_read_poscar(&mut stream).unwrap();

        let mut it = stream.lines().map(|l| l.unwrap());
        assert_eq!(it.next(), Some("    2    3    4".to_string()));
    }

    #[test]
    fn test_read_chg() {
        let mut stream = io::Cursor::new(SAMPLE.as_bytes());
        ChgBase::_read_poscar(&mut stream).unwrap();

        let chg = ChgBase::_read_chg(&mut stream).unwrap();
        assert_eq!(chg.shape(), &[2usize, 3, 4]);
        assert_eq!(chg[[1, 2, 3]], 0.10568153616E+01);
    }

}
