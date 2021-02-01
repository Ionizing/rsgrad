type MatX3<T> = Vec<[T;3]>;  // Nx3 matrix
type Mat33<T> = [[T;3];3];   // 3x3 matrix

use std::io;
use std::path::Path;
use std::fs;
use rayon::prelude;
use regex::Regex;


// DONE ISPIN
// DONE ions per type
// DONE element symbol
// DONE NKPTS
// DONE stress
// DONE cell
// DONE positions and forces
// DONE magmom
// DONE E-fermi
// DONE scf
// TODO viberation
// DONE LSORBIT
// DONE IBRION


#[derive(Clone)]
struct IonicIteration {
    pub nscf      : i32,
    pub toten     : f64,
    pub toten_z   : f64,
    pub cputime   : f64,
    pub magmom    : Option<Vec<f64>>,  // differs when ISPIN=1,2 and ncl versions
    pub positions : MatX3<f64>,
    pub forces    : MatX3<f64>,
    pub cell      : Mat33<f64>,
}

impl IonicIteration {
    pub fn new(nscf: i32, toten: f64, toten_z: f64, cputime: f64,
               magmom: Option<Vec<f64>>, positions: MatX3<f64>,
               forces: MatX3<f64>, cell: Mat33<f64>) -> Self {
        Self {
            nscf, toten, toten_z, cputime, magmom, positions, forces, cell
        }
    }
}


#[derive(Clone)]
struct Viberation {
    pub freq   : f64,  // in THz
    pub dxdydz : MatX3<f64>,
}

impl Viberation {
    pub fn new(freq: f64, dxdydz: MatX3<f64>) -> Self {
        Self {freq, dxdydz}
    }
}


#[derive(Clone)]
struct Outcar {
    pub lsorbit       : bool,
    pub ispin         : i32,
    pub ibrion        : i32,
    pub nions         : i32,
    pub nkpts         : i32,
    pub nbands        : i32,
    pub efermi        : f64,
    pub cell          : Mat33<f64>,
    pub ions_per_type : Vec<i32>,
    pub ion_types     : Vec<String>,
    pub scf           : Vec<IonicIteration>,
    pub vib           : Option<Viberation>,
}


impl Outcar {
    pub fn from_file(path: &(impl AsRef<Path> + ?Sized)) -> io::Result<Self> {
        let context: String = fs::read_to_string(path)?;

        let lsorbit         = Self::parse_lsorbit(&context);
        let ispin           = Self::parse_ispin(&context);
        let ibrion          = Self::parse_ibrion(&context);
        let nions           = Self::parse_nions(&context);
        let (nkpts, nbands) = Self::parse_nkpts_nbands(&context);
        let efermi          = Self::parse_efermi(&context);
        let cell            = Self::parse_cell(&context);
        let ions_per_type   = Self::parse_ions_per_type(&context);
        let ion_types       = Self::parse_ion_types(&context);

        let nscfv          = Self::parse_nscfs(&context);
        let totenv         = Self::parse_toten(&context);
        let toten_zv       = Self::parse_toten_z(&context);
        let cputimev       = Self::parse_cputime(&context);
        let (posv, forcev) = Self::parse_posforce(&context);
        let cellv          = Self::parse_opt_cells(&context);

        // Do some check
        let len = totenv.len();
        assert_eq!(nscfv.len()    , len);
        assert_eq!(toten_zv.len() , len);
        assert_eq!(cputimev.len() , len);
        assert_eq!(posv.len()     , len);
        assert_eq!(forcev.len()   , len);
        assert_eq!(cellv.len()    , len);

        // let scf = (0..len).fold(vec![], |mut iscf, i| {
        //     let ioniter = IonicIteration::new(nscfv[i], totenv[i], toten_zv[i], cputimev[i],
        //                                       None, posv[i].clone(), forcev[i].clone(), cellv[i]);
        //     iscf.push(ioniter);
        //     iscf
        // });

        let scf = nscfv.into_iter()
                       .zip(totenv.into_iter())
                       .zip(toten_zv.into_iter())
                       .zip(cputimev.into_iter())
                       .zip(posv.into_iter())
                       .zip(forcev.into_iter())
                       .zip(cellv.into_iter())
                       .map(|((((((iscf, e), ez), cpu), pos), f), cell)| {
                           IonicIteration::new(iscf, e, ez, cpu, None, pos, f, cell)
                       })
                       .collect::<Vec<IonicIteration>>();

        let vib = None;

        Ok(
            Outcar {
                lsorbit,
                ispin,
                ibrion,
                nions,
                nkpts,
                nbands,
                efermi,
                cell,
                ions_per_type,
                ion_types,
                scf,
                vib
            }
        )
    }

    fn parse_ispin(context: &str) -> i32 {
        Regex::new(r"ISPIN  =      (\d)")
            .unwrap()
            .captures(context)
            .unwrap()
            .get(1)
            .unwrap()
            .as_str()
            .parse::<i32>()
            .unwrap()
    }

    fn parse_nions(context: &str) -> i32 {
        Regex::new(r"NIONS = \s+(\d+)")
            .unwrap()
            .captures(context)
            .unwrap()
            .get(1)
            .unwrap()
            .as_str()
            .parse::<i32>()
            .unwrap()
    }

    fn parse_toten(context: &str) -> Vec<f64> {
        Regex::new(r"free  energy   TOTEN  = \s*([-+]?[0-9]+[.]?[0-9]*([eE][-+]?[0-9]+)?) eV")
            .unwrap()
            .captures_iter(context)
            .map(|x| {
                x.get(1)
                 .unwrap()
                 .as_str()
                 .parse::<f64>()
                    .unwrap()
            })
            .collect()
    }

    fn parse_toten_z(context: &str) -> Vec<f64> {
        Regex::new(r"energy  without entropy=\s+([-+]?[0-9]+[.]?[0-9]*?)  energy\(sigma->0\) =\s+([-+]?[0-9]+[.]?[0-9]*)")
            .unwrap()
            .captures_iter(context)
            // .inspect(|x| {dbg!(&x);})
            .map(|x| {
                x.get(2)
                 .unwrap()
                 .as_str()
                 .parse::<f64>()
                    .unwrap()
            })
            .collect()
    }

    fn parse_cputime(context: &str) -> Vec<f64> {
        Regex::new(r"LOOP\+:  cpu time .* real time \s+([0-9]+[.]?[0-9]*)")
            .unwrap()
            .captures_iter(context)
            .map(|x| {
                x.get(1)
                 .unwrap()
                 .as_str()
                 .parse::<f64>()
                    .unwrap()
            })
            .collect()
    }

    fn parse_magmoms(context: &str) -> Option<Vec<Vec<f64>>> {
        Regex::new(r"free  energy")
            .unwrap()
            .find_iter(context)
            .map(|x| x.start())
            .map(|x| Self::_parse_magmom(&context[..x]))
            .inspect(|x| {dbg!(x);})
            .collect()
    }

    fn _parse_magmom(context: &str) -> Option<Vec<f64>> {
        let pos = context
            .rmatch_indices("number of electron")
            .next()
            .unwrap()
            .0;
        let ret = context[pos..]
            .lines()
            .next()
            .unwrap()
            .split_whitespace()
            .skip(5)
            .map(|x| x.trim().parse::<f64>().unwrap())
            .collect::<Vec<_>>();
        match ret.len() {
            0 => None,
            _ => Some(ret)
        }
    }

    fn parse_posforce(context: &str) -> (Vec<MatX3<f64>>, Vec<MatX3<f64>>) {
        Regex::new(r"(?m)^ POSITION")
            .unwrap()
            .find_iter(context)
            .map(|x| x.start())
            .map(|x| {
                Self::_parse_posforce_single_iteration(&context[x..])
            })
            .fold((vec![], vec![]), |mut acc, (p, f)| {
                acc.0.push(p);
                acc.1.push(f);
                acc
            })
    }

    fn _parse_posforce_single_iteration(context: &str) -> (MatX3<f64>, MatX3<f64>) {
        assert!(context.starts_with(" POSITION"));
        context.lines()
               // .inspect(|x| {dbg!(&x);})
               .skip(2)
               .take_while(|x| !x.starts_with(" ----"))
               // .inspect(|x| {dbg!(&x);})
               .map(|x| {
                   x.split_whitespace()
                       .map(|x| x.parse::<f64>().unwrap())
                       .collect::<Vec<f64>>()
               })
               .fold((vec![], vec![]), |mut ret, x|{
                   ret.0.push([x[0], x[1], x[2]]);
                   ret.1.push([x[3], x[4], x[5]]);
                   ret
               })
    }

    fn parse_efermi(context: &str) -> f64 {
        Regex::new(r" E-fermi : \s+([-+]?[0-9]+[.]?[0-9]*)")
            .unwrap()
            .captures(context)
            .unwrap()
            .get(1)
            .unwrap()
            .as_str()
            .parse::<f64>()
            .unwrap()
    }

    fn parse_nkpts_nbands(context: &str) -> (i32, i32) {
        let v = Regex::new(r"NKPTS = \s*(\d+) .* NBANDS= \s*(\d+)")
            .unwrap()
            .captures(context)
            .unwrap()
            .iter()
            .skip(1)
            .map(|x| {
                x.unwrap()
                 .as_str()
                 .parse::<i32>()
                    .unwrap()
            })
            .collect::<Vec<i32>>();
        (v[0], v[1])
    }

    fn parse_cell(context: &str) -> Mat33<f64> {
        let pos = Regex::new(r"direct lattice vectors")
            .unwrap()
            .find(context)
            .unwrap()
            .start();
        let v = &context[pos..]
            .lines()
            .skip(1)
            .take(3)
            .map(|l| {
                let v = l.split_whitespace()
                         .map(|x| x.parse::<f64>().unwrap())
                         .collect::<Vec<f64>>();
                [v[0], v[1], v[2]]
            })
            .collect::<Vec<[f64; 3]>>();
        [v[0], v[1], v[2]]
    }

    fn parse_opt_cells(context: &str) -> Vec<Mat33<f64>> {
        Regex::new(r"direct lattice vectors")
            .unwrap()
            .find_iter(context)
            .map(|x| x.start())
            .skip(1)
            .map(|x| Self::parse_cell(&context[x..]))
            .collect()
    }

    fn parse_ions_per_type(context: &str) -> Vec<i32> {
        Regex::new(r"(?m)ions per type = .*$")
            .unwrap()
            .find(context)
            .unwrap()
            .as_str()
            .split_whitespace()
            .skip(4)
            .map(|x| x.parse::<i32>().unwrap())
            .collect()
    }

    fn parse_ion_types(context: &str) -> Vec<String> {
        let mut v = Regex::new(r"(?m)^ POTCAR:.*$")
            .unwrap()
            .find_iter(context)
            .map(|l| {
                l.as_str()
                 .split_whitespace()
                 .nth(2)
                 .unwrap()
                 .to_owned()
            })
            .collect::<Vec<String>>();

        let len = v.len() / 2;
        (0..len).for_each(|_| {v.pop();});
        v
    }

    fn parse_nscfs(context: &str) -> Vec<i32> {
        Regex::new(r"free  energy")  // navigate to tail of ionic step
            .unwrap()
            .find_iter(context)
            .map(|x| x.start())
            .map(|x| Self::_parse_nscf(&context[..x]))
            .collect()
    }

    fn _parse_nscf(context: &str) -> i32 {
        let pos = context
            .rmatch_indices("Iteration") // get the last "Iteration" during ionic step
            .next()
            .unwrap()
            .0;
        let context = &context[pos..];
        Regex::new(r"Iteration\s*\d+\(\s*(\d+)\)")
            .unwrap()
            .captures(context)
            .unwrap()
            .get(1)
            .unwrap()
            .as_str()
            .parse::<i32>()
            .unwrap()
    }

    fn parse_stress(context: &str) -> Vec<f64> {
        Regex::new(r"external pressure = \s*(\S+) kB")
            .unwrap()
            .captures_iter(context)
            .map(|x| {
                x.get(1)
                 .unwrap()
                 .as_str()
                 .parse::<f64>()
                    .unwrap()
            })
            .collect()
    }

    fn parse_ibrion(context: &str) -> i32 {
        Regex::new(r"IBRION = \s*(\S+) ")
            .unwrap()
            .captures(context)
            .unwrap()
            .get(1)
            .unwrap()
            .as_str()
            .parse::<i32>()
            .unwrap()
    }

    fn parse_lsorbit(context: &str) -> bool {
        match Regex::new(r"LSORBIT\s*=\s*([TF])")
            .unwrap()
            .captures(context)
            .unwrap()
            .get(1)
            .unwrap()
            .as_str() {
                "T" => true,
                "F" => false,
                _ => unreachable!("Invalid value for LSORBIT, should be T or F")
            }
    }
}


#[cfg(test)]
mod tests{
    use super::*;

    #[test]
    fn test_parse_ispin() {
        let input = r#"
   ICHARG =      2    charge: 1-file 2-atom 10-const
   ISPIN  =      1    spin polarized calculation?
   LNONCOLLINEAR =      F non collinear calculations"#;
        assert_eq!(Outcar::parse_ispin(&input), 1i32);
    }

    #[test]
    fn test_parse_nions() {
        let input = r#"
   k-points           NKPTS =      1   k-points in BZ     NKDIM =      1   number of bands    NBANDS=      8
   number of dos      NEDOS =    301   number of ions     NIONS =      4
   non local maximal  LDIM  =      4   non local SUM 2l+1 LMDIM =      8 "#;
        assert_eq!(Outcar::parse_nions(&input), 4i32);
    }

    #[test]
    fn test_parse_toten() {
        let input = r#"
  free energy    TOTEN  =        51.95003235 eV
  free energy    TOTEN  =       -10.91478741 eV
  free energy    TOTEN  =       -22.11911831 eV
  free  energy   TOTEN  =       -19.26550806 eV
  free  energy   TOTEN  =       -19.25519593 eV
  free  energy   TOTEN  =       -19.26817124 eV
"#;
        let output = vec![-19.26550806f64, -19.25519593, -19.26817124];
        assert_eq!(Outcar::parse_toten(&input), output);
    }

    #[test]
    fn test_parse_toten_z() {
        let input = r#"
  energy without entropy =       51.93837380  energy(sigma->0) =       51.94614617
  energy without entropy =      -10.92638322  energy(sigma->0) =      -10.91865268
  energy without entropy =      -22.13071412  energy(sigma->0) =      -22.12298358
  energy  without entropy=      -19.27710387  energy(sigma->0) =      -19.26937333
  energy  without entropy=      -19.26679174  energy(sigma->0) =      -19.25906120
  energy  without entropy=      -19.27976705  energy(sigma->0) =      -19.27203651"#;
        let output = vec![-19.26937333f64, -19.25906120, -19.27203651];
        assert_eq!(Outcar::parse_toten_z(&input), output);
    }

    #[test]
    fn test_parse_cputime() {
        let input = r#"
      LOOP:  cpu time    0.0894: real time    0.0949
      LOOP:  cpu time    0.0360: real time    0.0330
      LOOP:  cpu time    0.0275: real time    0.0261
     LOOP+:  cpu time    2.0921: real time    2.0863
     LOOP+:  cpu time    1.2021: real time    1.1865
     LOOP+:  cpu time    1.2788: real time    1.2670"#;
        let output = vec![2.0863, 1.1865, 1.2670];
        assert_eq!(Outcar::parse_cputime(&input), output);
    }

    #[test]
    fn test_parse_posforce_single_iteration() {
        let input = r#" POSITION                                       TOTAL-FORCE (eV/Angst)
 -----------------------------------------------------------------------------------
      3.87720      4.01520      4.00000        -0.438233     -0.328151      0.000000
      3.00000      2.48290      4.00000         0.000000      0.536218      0.000000
      2.12280      4.01520      4.00000         0.438233     -0.328151      0.000000
      3.00000      3.50000      4.00000         0.000000      0.120085      0.000000
 -----------------------------------------------------------------------------------
    total drift:                                0.000000     -0.000260     -0.000000 "#;
        let output = (
            vec![[3.87720,4.01520,4.00000],
                 [3.00000,2.48290,4.00000],
                 [2.12280,4.01520,4.00000],
                 [3.00000,3.50000,4.00000]]
            ,
            vec![[-0.438233, -0.328151, 0.000000],
                 [ 0.000000,  0.536218, 0.000000],
                 [ 0.438233, -0.328151, 0.000000],
                 [ 0.000000,  0.120085, 0.000000]]
        );

        assert_eq!(Outcar::_parse_posforce_single_iteration(&input), output);
    }

    #[test]
    fn test_parse_posforce() {
        let input = r#"
 POSITION                                       TOTAL-FORCE (eV/Angst)
 -----------------------------------------------------------------------------------
      3.87720      4.01520      4.00000        -0.438233     -0.328151      0.000000
      3.00000      2.48290      4.00000         0.000000      0.536218      0.000000
      2.12280      4.01520      4.00000         0.438233     -0.328151      0.000000
      3.00000      3.50000      4.00000         0.000000      0.120085      0.000000
 -----------------------------------------------------------------------------------
--
 POSITION                                       TOTAL-FORCE (eV/Angst)
 -----------------------------------------------------------------------------------
      3.89220      4.01520      4.00000        -0.930834     -0.563415      0.000000
      3.00000      2.48290      4.00000        -0.006828      0.527001      0.000000
      2.12280      4.01520      4.00000         0.458533     -0.304111      0.000000
      3.00000      3.50000      4.00000         0.479129      0.340525      0.000000
 -----------------------------------------------------------------------------------
--
 POSITION                                       TOTAL-FORCE (eV/Angst)
 -----------------------------------------------------------------------------------
      3.86220      4.01520      4.00000         0.089245     -0.065055      0.000000
      3.00000      2.48290      4.00000         0.007618      0.545925      0.000000
      2.12280      4.01520      4.00000         0.417195     -0.352508      0.000000
      3.00000      3.50000      4.00000        -0.514057     -0.128362      0.000000
 -----------------------------------------------------------------------------------
--"#;
        let output = (
            vec![vec![[3.87720,4.01520,4.00000],
                      [3.00000,2.48290,4.00000],
                      [2.12280,4.01520,4.00000],
                      [3.00000,3.50000,4.00000]],
                 vec![[3.89220,4.01520,4.00000],
                      [3.00000,2.48290,4.00000],
                      [2.12280,4.01520,4.00000],
                      [3.00000,3.50000,4.00000]],
                 vec![[3.86220,4.01520,4.00000],
                      [3.00000,2.48290,4.00000],
                      [2.12280,4.01520,4.00000],
                      [3.00000,3.50000,4.00000]],
            ],

            vec![vec![[-0.438233, -0.328151, 0.000000],
                      [ 0.000000,  0.536218, 0.000000],
                      [ 0.438233, -0.328151, 0.000000],
                      [ 0.000000,  0.120085, 0.000000]],
                 vec![[-0.930834, -0.563415, 0.000000],
                      [-0.006828,  0.527001, 0.000000],
                      [ 0.458533, -0.304111, 0.000000],
                      [ 0.479129,  0.340525, 0.000000]],
                 vec![[ 0.089245, -0.065055, 0.000000],
                      [ 0.007618,  0.545925, 0.000000],
                      [ 0.417195, -0.352508, 0.000000],
                      [-0.514057, -0.128362, 0.000000]]
            ]
        );
        assert_eq!(Outcar::parse_posforce(&input), output);
    }

    #[test]
    fn test_parse_efermi() {
        let input = " E-fermi :  -0.7865     XC(G=0):  -2.0223     alpha+bet : -0.5051";
        let output = -0.7865f64;
        assert_eq!(Outcar::parse_efermi(&input), output);
    }

    #[test]
    fn test_parse_nkpts_nbands() {
        let input = r#"
 Dimension of arrays:
   k-points           NKPTS =      1   k-points in BZ     NKDIM =      1   number of bands    NBANDS=      8
   number of dos      NEDOS =    301   number of ions     NIONS =      4"#;
        let output = (1i32, 8i32);
        assert_eq!(Outcar::parse_nkpts_nbands(&input), output);
    }

    #[test]
    fn test_parse_cell() {
        let input = r#"
  energy-cutoff  :      400.00
  volume of cell :      336.00
      direct lattice vectors                 reciprocal lattice vectors
     6.000000000  0.000000000  0.000000000     0.166666667  0.000000000  0.000000000
     0.000000000  7.000000000  0.000000000     0.000000000  0.142857143  0.000000000
     0.000000000  0.000000000  8.000000000     0.000000000  0.000000000  0.125000000 "#;
        let output = [[6.0, 0.0, 0.0],
                      [0.0, 7.0, 0.0],
                      [0.0, 0.0, 8.0]];
        assert_eq!(Outcar::parse_cell(&input), output);
    }

    #[test]
    fn test_parse_opt_cells() {
        let input = r#"
      direct lattice vectors                 reciprocal lattice vectors
     6.000000000  0.000000000  0.000000000     0.166666667  0.000000000  0.000000000
     0.000000000  7.200000000  0.000000000     0.000000000  0.142857143  0.000000000
     0.000000000  0.000000000  8.000000000     0.000000000  0.000000000  0.125000000
--
      direct lattice vectors                 reciprocal lattice vectors
     6.000000000  0.000000000  0.000000000     0.166666667  0.000000000  0.000000000
     0.000000000  7.000000000  0.000000000     0.000000000  0.142857143  0.000000000
     0.000000000  0.000000000  8.000000000     0.000000000  0.000000000  0.125000000
--
      direct lattice vectors                 reciprocal lattice vectors
     6.000000000  0.000000000  0.000000000     0.166666667  0.000000000  0.000000000
     0.000000000  7.000000000  0.000000000     0.000000000  0.142857143  0.000000000
     0.000000000  0.000000000  8.000000000     0.000000000  0.000000000  0.125000000
--"#;
        let output = vec![ [[6.0, 0.0, 0.0],
                            [0.0, 7.0, 0.0],
                            [0.0, 0.0, 8.0]]; 2];
        assert_eq!(Outcar::parse_opt_cells(&input), output);
    }

    #[test]
    fn test_parse_ions_per_type() {
        let input = r#"
   support grid    NGXF=    60 NGYF=   72 NGZF=   80
   ions per type =               3   1
 NGX,Y,Z   is equivalent  to a cutoff of   8.31,  8.55,  8.31 a.u. "#;
        let output = vec![3i32, 1];
        assert_eq!(Outcar::parse_ions_per_type(&input), output);
    }


    #[test]
    fn test_parse_ion_types() {
        let input = r#"
 INCAR:
 POTCAR:    PAW_PBE H 15Jun2001
 POTCAR:    PAW_PBE N 08Apr2002
 POTCAR:    PAW_PBE H 15Jun2001
   VRHFIN =H: ultrasoft test
   LEXCH  = PE
   EATOM  =    12.4884 eV,    0.9179 Ry
......
   number of l-projection  operators is LMAX  =           3
   number of lm-projection operators is LMMAX =           5

 POTCAR:    PAW_PBE N 08Apr2002
   VRHFIN =N: s2p3
   LEXCH  = PE
   EATOM  =   264.5486 eV,   19.4438 Ry"#;
        let output = vec!["H", "N"];
        assert_eq!(Outcar::parse_ion_types(&input), output);
    }


    #[test]
    fn test_parse_nscf() {
        let input = r#"
 -0.508   0.103   0.035   0.000   0.070
----------------------------------------- Iteration    1(  23)  ---------------------------------------
......

  FREE ENERGIE OF THE ION-ELECTRON SYSTEM (eV)
  ---------------------------------------------------
  free  energy   TOTEN  =       -19.26550806 eV
  energy  without entropy=      -19.27710387  energy(sigma->0) =      -19.26937333 "#;
        let output = 23i32;
        assert_eq!(Outcar::_parse_nscf(&input), output);
    }

    #[test]
    fn test_parse_nscfs() {
        let input = r#"
----------------------------------------- Iteration    1(  22)  ---------------------------------------
----------------------------------------- Iteration    1(  23)  ---------------------------------------
......

  FREE ENERGIE OF THE ION-ELECTRON SYSTEM (eV)
  ---------------------------------------------------
  free  energy   TOTEN  =       -19.26550806 eV
  energy  without entropy=      -19.27710387  energy(sigma->0) =      -19.26937333
......


----------------------------------------- Iteration    2(  12)  ---------------------------------------
----------------------------------------- Iteration    2(  13)  ---------------------------------------
......
  FREE ENERGIE OF THE ION-ELECTRON SYSTEM (eV)
  ---------------------------------------------------
  free  energy   TOTEN  =       -19.25519593 eV
  energy  without entropy=      -19.26679174  energy(sigma->0) =      -19.25906120
......


----------------------------------------- Iteration    3(  12)  ---------------------------------------
----------------------------------------- Iteration    3(  13)  ---------------------------------------
......
  FREE ENERGIE OF THE ION-ELECTRON SYSTEM (eV)
  ---------------------------------------------------
  free  energy   TOTEN  =       -19.26817124 eV
  energy  without entropy=      -19.27976705  energy(sigma->0) =      -19.27203651
"#;
        let output = vec![23, 13, 13];
        assert_eq!(Outcar::parse_nscfs(&input), output);
    }

    #[test]
    fn test_parse_stress() {
        let input = r#"
  in kB      -6.78636    -7.69902    -4.03340     0.00000     0.00000     0.00000
  external pressure =       -6.17 kB  Pullay stress =        0.00 kB
--
  in kB      -8.92250    -8.14636    -4.01885    -1.10430     0.00000     0.00000
  external pressure =       -7.03 kB  Pullay stress =        0.00 kB
--
  in kB      -4.56989    -7.18734    -4.04843     1.18589     0.00000     0.00000
  external pressure =       -5.27 kB  Pullay stress =        0.00 kB"#;
        let output = vec![-6.17, -7.03, -5.27];
        assert_eq!(Outcar::parse_stress(&input), output);
    }

    #[test]
    fn test_parse_ibrion() {
        let input = r#"
   NSW    =     85    number of steps for IOM
   NBLOCK =      1;   KBLOCK =      1    inner block; outer block
   IBRION =      5    ionic relax: 0-MD 1-quasi-New 2-CG
   NFREE  =      2    steps in history (QN), initial steepest desc. (CG)
   ISIF   =      2    stress and relaxation
"#;
        let output = 5i32;
        assert_eq!(Outcar::parse_ibrion(&input), output);
    }

    #[test]
    fn test_parse_lsorbit() {
        let input = r#"
   LNONCOLLINEAR =      F non collinear calculations
   LSORBIT =      F    spin-orbit coupling
   INIWAV =      1    electr: 0-lowe 1-rand  2-diag "#;
        let output = false;
        assert_eq!(Outcar::parse_lsorbit(&input), output);
    }


    #[test]
    fn test_parse_magmoms() {
        let input = r#"
 total energy-change (2. order) :-0.5897058E-05  (-0.8072299E-08)
 number of electron     309.9999998 magnetization      42.0005098
 augmentation part       88.5937960 magnetization      26.8073410
......
  free  energy   TOTEN  =      -391.79003630 eV
  energy  without entropy=     -391.77828290  energy(sigma->0) =     -391.78611850
"#;
        let output = Some(vec![vec![42.0005098f64]]);
        assert_eq!(Outcar::parse_magmoms(&input), output);


        let input = r#"
 total energy-change (2. order) :-0.5897058E-05  (-0.8072299E-08)
 number of electron     309.9999998 magnetization      42.0005098 42.0005098 42.0005098
 augmentation part       88.5937960 magnetization      26.8073410 26.8073410 26.8073410
......
  free  energy   TOTEN  =      -391.79003630 eV
  energy  without entropy=     -391.77828290  energy(sigma->0) =     -391.78611850
"#;
        let output = Some(vec![vec![42.0005098f64; 3]]);
        assert_eq!(Outcar::parse_magmoms(&input), output);


        let input = r#"
 total energy-change (2. order) :-0.5897058E-05  (-0.8072299E-08)
 number of electron     309.9999998 magnetization
 augmentation part       88.5937960 magnetization
......
  free  energy   TOTEN  =      -391.79003630 eV
  energy  without entropy=     -391.77828290  energy(sigma->0) =     -391.78611850
"#;
        let output = None;
        assert_eq!(Outcar::parse_magmoms(&input), output);


        let input = r#"
 total energy-change (2. order) :-0.5897058E-05  (-0.8072299E-08)
 number of electron     309.9999998 magnetization      42.0005098 42.0005098 42.0005098
 augmentation part       88.5937960 magnetization      26.8073410 26.8073410 26.8073410
......
  free  energy   TOTEN  =      -391.79003630 eV
  energy  without entropy=     -391.77828290  energy(sigma->0) =     -391.78611850

 total energy-change (2. order) :-0.5897058E-05  (-0.8072299E-08)
 number of electron     309.9999998 magnetization      42.0005098 42.0005098 42.0005098
 augmentation part       88.5937960 magnetization      26.8073410 26.8073410 26.8073410
......
  free  energy   TOTEN  =      -391.79003630 eV
  energy  without entropy=     -391.77828290  energy(sigma->0) =     -391.78611850

 total energy-change (2. order) :-0.5897058E-05  (-0.8072299E-08)
 number of electron     309.9999998 magnetization      42.0005098 42.0005098 42.0005098
 augmentation part       88.5937960 magnetization      26.8073410 26.8073410 26.8073410
......
  free  energy   TOTEN  =      -391.79003630 eV
  energy  without entropy=     -391.77828290  energy(sigma->0) =     -391.78611850
"#;
        let output = Some(vec![vec![42.0005098f64; 3]; 3]);
        assert_eq!(Outcar::parse_magmoms(&input), output);
    }
}
