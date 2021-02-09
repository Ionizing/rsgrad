type MatX3<T> = Vec<[T;3]>;  // Nx3 matrix
type Mat33<T> = [[T;3];3];   // 3x3 matrix

use std::io;
use std::path::Path;
use std::fs;
use std::fmt;
use rayon;
use regex::Regex;
use itertools::multizip;
use colored::Colorize;

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
// DONE viberation
// DONE LSORBIT
// DONE IBRION
// DONE ion masses


#[derive(Clone, PartialEq, Debug)]
pub struct IonicIteration {
    pub nscf      : i32,
    pub toten     : f64,
    pub toten_z   : f64,
    pub cputime   : f64,
    pub stress    : f64,
    pub magmom    : Option<Vec<f64>>,  // differs when ISPIN=1,2 and ncl versions
    pub positions : MatX3<f64>,
    pub forces    : MatX3<f64>,
    pub cell      : Mat33<f64>,
}

impl IonicIteration {
    pub fn new(nscf: i32, toten: f64, toten_z: f64, cputime: f64,
               stress: f64, magmom: Option<Vec<f64>>, positions: MatX3<f64>,
               forces: MatX3<f64>, cell: Mat33<f64>) -> Self {
        Self {
            nscf, toten, toten_z, cputime, stress,
            magmom, positions, forces, cell
        }
    }
    // The parsing process is done within `impl Outcar`
}


impl fmt::Display for IonicIteration {
    fn fmt(&self, f: &mut fmt::Formatter<'_>) -> fmt::Result {
        let nscf_text    = format!("{:4}", self.nscf).normal();
        let toten_text   = format!("{:11.5}", self.toten).normal();
        let totez_text   = format!("{:11.5}", self.toten_z).bright_green();
        let cputime_text = format!("{:6.2}", self.cputime / 60.0).white();
        let stress_text  = format!("{:6.2}", self.stress).normal();

        let fsize = self.forces.iter()
                               .map(|f| (f[0]*f[0] + f[1]*f[1] * f[2]*f[2]).sqrt())
                               .collect::<Vec<_>>();
        let maxf_text   = format!("{:6.3}", fsize.iter().cloned().fold(0.0, f64::max)).bright_green();
        let avgf_text   = format!("{:6.3}", fsize.iter().sum::<f64>() / self.forces.len() as f64).normal();
        let magmom_text = if let Some(mag) = &self.magmom {
            mag.iter()
               .map(|n| format!("{:8.4}", n))
               .collect::<Vec<_>>()
                .join(" ")
                .bright_yellow()
        } else { "   NoMag".to_string().normal() };

        write!(f, "{E} {Ez} {nscf} {Favg} {Fmax} {Stress} {Time} {Mag}",
               E=toten_text, Ez=totez_text, nscf=nscf_text, Favg=avgf_text,
               Fmax=maxf_text, Time=cputime_text, Mag=magmom_text, Stress=stress_text)
    }
}


pub struct PrintOptIterations(Vec<IonicIteration>);

impl fmt::Display for PrintOptIterations {
    fn fmt(&self, f: &mut fmt::Formatter<'_>) -> fmt::Result {
        let header = format!("# NStep {:>11} {:>11} {:>4} {:>6} {:>6} {:>6} {:>6} {:>8}",
               "E/eV", "Ez/eV", "NSCF", "Favg", "Fmax", "Stress", "t/min", "Mag/muB");
        let body = self.0.iter().enumerate()
            .map(|(i, x)| format!(" {:6} {}", i+1, x))
            .collect::<Vec<String>>()
            .join("\n");

        write!(f, "{}\n{}", header.bright_green(), body)
    }
}

impl From<Vec<IonicIteration>> for PrintOptIterations {
    fn from(v: Vec<IonicIteration>) -> Self {
        Self(v)
    }
}


#[derive(Clone, PartialEq, Debug)]
pub struct Viberation {
    pub freq       : f64,  // in THz
    pub dxdydz     : MatX3<f64>,
    pub is_imagine : bool, // denote wheher this mode is an imagine mode
}

impl Viberation {
    pub fn new(freq: f64, dxdydz: MatX3<f64>, is_imagine: bool) -> Self {
        Self {freq, dxdydz, is_imagine}
    }
    // The parsing process is done withon `impl Outcar`
}


#[derive(Clone, Debug, PartialEq)]
pub struct Outcar {
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
    pub ion_masses    : Vec<f64>,  // .len() == nions
    pub ion_iters     : Vec<IonicIteration>,
    pub vib           : Option<Vec<Viberation>>, // .len() == degrees of freedom
}


impl Outcar {
    pub fn from_file(path: &(impl AsRef<Path> + ?Sized)) -> io::Result<Self> {
        let context: String = fs::read_to_string(path)?;

        let mut lsorbit         = false;
        let mut ispin           = 0i32;
        let mut ibrion          = 0i32;
        let mut nions           = 0i32;
        let (mut nkpts, mut nbands) = (0i32, 0i32);
        let mut efermi          = 0.0f64;
        let mut cell            = [[0.0f64; 3]; 3];
        let mut ext_pressure    = vec![0.0f64; 0];
        let mut ions_per_type   = vec![0i32; 0];
        let mut ion_types       = vec![String::new();0];
        let mut ion_masses      = vec![0.0f64; 0];

        let mut nscfv          = vec![0i32; 0];
        let mut totenv         = vec![0.0f64; 0];
        let mut toten_zv       = vec![0.0f64; 0];
        let mut magmomv        = vec![Some(vec![0.0f64; 0]); 0];
        let mut cputimev       = vec![0.0f64; 0];
        let (mut posv, mut forcev) = (vec![vec![[0.0f64; 3];0]; 0], vec![vec![[0.0f64; 3];0]; 0]);
        let mut cellv          = vec![[[0.0f64; 3]; 3]; 0];

        rayon::scope(|s| {
            s.spawn(|_| { lsorbit         = Self::parse_lsorbit(&context) });
            s.spawn(|_| { ispin           = Self::parse_ispin(&context) });
            s.spawn(|_| { ibrion          = Self::parse_ibrion(&context) });
            s.spawn(|_| { nions           = Self::parse_nions(&context) });
            s.spawn(|_| {
                let (_nkpts, _nbands) = Self::parse_nkpts_nbands(&context);
                nkpts = _nkpts;
                nbands = _nbands;
            });
            s.spawn(|_| { efermi          = Self::parse_efermi(&context) });
            s.spawn(|_| { cell            = Self::parse_cell(&context) });
            s.spawn(|_| { ext_pressure    = Self::parse_stress(&context) });
            s.spawn(|_| { ions_per_type   = Self::parse_ions_per_type(&context) });
            s.spawn(|_| { ion_types       = Self::parse_ion_types(&context) });
            s.spawn(|_| { ion_masses      = Self::parse_ion_masses(&context) });

            s.spawn(|_| { nscfv          = Self::parse_nscfs(&context) });
            s.spawn(|_| { totenv         = Self::parse_toten(&context) });
            s.spawn(|_| { toten_zv       = Self::parse_toten_z(&context) });
            s.spawn(|_| { magmomv        = Self::parse_magmoms(&context) });
            s.spawn(|_| { cputimev       = Self::parse_cputime(&context) });
            s.spawn(|_| {
                let (_posv, _forcev) = Self::parse_posforce(&context);
                posv = _posv;
                forcev = _forcev;
            });
            s.spawn(|_| { cellv          = Self::parse_opt_cells(&context) });
        });

        // Do some check
        let len = totenv.len();
        assert_eq!(nscfv.len()    , len);
        assert_eq!(toten_zv.len() , len);
        assert_eq!(cputimev.len() , len);
        assert_eq!(posv.len()     , len);
        assert_eq!(forcev.len()   , len);
        assert_eq!(cellv.len()    , len);

        let ion_iters = multizip((nscfv, totenv, toten_zv, magmomv, cputimev, ext_pressure, posv, forcev, cellv))
            .map(|(iscf, e, ez, mag, cpu, stress, pos, f, cell)| {
                IonicIteration::new(iscf, e, ez, cpu, stress, mag, pos, f, cell)
            })
            .collect::<Vec<IonicIteration>>();

        let vib = Self::parse_viberations(&context);

        Ok(
            Self {
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
                ion_masses,
                ion_iters,
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
        Regex::new(r"free  energy   TOTEN  = \s*(\S+) eV")
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
        Regex::new(r"energy  without entropy=\s+(?:\S+)  energy\(sigma->0\) =\s+(\S+)")
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

    fn parse_cputime(context: &str) -> Vec<f64> {
        Regex::new(r"LOOP\+:  cpu time .* real time\s*(\S+)")
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

    fn parse_magmoms(context: &str) -> Vec<Option<Vec<f64>>> {
        Regex::new(r"free  energy")
            .unwrap()
            .find_iter(context)
            .map(|x| x.start())
            .map(|x| Self::_parse_magmom(&context[..x]))
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
        Regex::new(r"(?m)^ POSITION \s+ TOTAL-FORCE \(eV/Angst\)")
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
               .skip(2)
               .take_while(|x| !x.starts_with(" ----"))
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
        Regex::new(r" E-fermi : \s+(\S+)")
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
        let skip_cnt: usize = if context.find(" old parameters").is_some() {
            2
        } else {
            1
        };
        Regex::new(r"direct lattice vectors")
            .unwrap()
            .find_iter(context)
            .map(|x| x.start())
            .skip(skip_cnt)
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

    fn parse_ion_masses(context: &str) -> Vec<f64> {
        let ions_per_type = Self::parse_ions_per_type(context);
        let masses_per_type = Regex::new(r"POMASS = \s*(\S+); ZVAL")
            .unwrap()
            .captures_iter(context)
            .map(|x| { x.get(1)
                       .unwrap()
                       .as_str()
                       .parse::<f64>()
                       .unwrap()
            })
            .collect::<Vec<f64>>();

        ions_per_type.into_iter()
            .zip(masses_per_type.into_iter())
            .fold(vec![], |mut acc, (n, m): (i32, f64)| {
                (0..n).for_each(|_| acc.push(m));
                acc
            })
    }

    fn parse_viberations(context: &str) -> Option<Vec<Viberation>> {
        let massess_sqrt = Self::parse_ion_masses(context)
            .iter()
            .map(|x| x.sqrt())
            .collect::<Vec<_>>();

        let ndof = Self::_parse_dof(context)? as usize;

        let mut vibs = Regex::new(r"(?m) .* 2PiTHz.* cm-1")
            .unwrap()
            .find_iter(context)
            .take(ndof)
            .map(|x| x.start())
            .map(|x| Self::_parse_single_vibmode(&context[x..]))
            .collect::<Vec<_>>();

        if vibs.is_empty() { return None; }

        vibs.iter_mut()
            .for_each(|v| {
                v.dxdydz.iter_mut()
                        .zip(massess_sqrt.iter())
                        .for_each(|(d, m)| {
                            d.iter_mut()
                             .for_each(|x| *x /= m)
                        })
            });

        Some(vibs)
    }

    fn _parse_single_vibmode(context: &str) -> Viberation {
        let freq = Regex::new(r"2PiTHz \s*(\S*) cm-1")
            .unwrap()
            .captures(context)
            .unwrap()
            .get(1)
            .unwrap()
            .as_str()
            .parse::<f64>()
            .unwrap();

        let is_imagine = match Regex::new(r"f(/i|  )= .* THz")  // Find the line contains "f/i=  xxxx THz"
            .unwrap()
            .captures(context)
            .unwrap()
            .get(1)
            .unwrap()
            .as_str() {
                "  " => false,
                "/i" => true,
                _ => unreachable!("Invalid viberation frequency indicator")
            };


        let start_pos = Regex::new(r"dx \s* dy \s* dz")
            .unwrap()
            .find(context)
            .unwrap()
            .start();

        let dxdydz: MatX3<f64> = context[start_pos..]
            .lines()
            .skip(1)
            .take_while(|l| !l.trim().is_empty())
            .map(|l| {
                let v = l.split_whitespace()
                         .skip(3)
                         .take(3)
                         .map(|token| token.parse::<f64>().unwrap())
                         .collect::<Vec<_>>();
                [v[0], v[1], v[2]]
            })
            .collect::<MatX3<f64>>();

        Viberation::new(freq, dxdydz, is_imagine)
    }

    fn _parse_dof(context: &str) -> Option<i32> {
        Regex::new(r"(?m)^   Degrees of freedom DOF   = \s*(\S+)$")
            .unwrap()
            .captures(context)?
            .get(1)
            .unwrap()
            .as_str()
            .parse::<i32>()
            .ok()
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
     LOOP+:  cpu time 1543.2679: real time 1544.6603
     LOOP+:  cpu time    1.2788: real time    1.2670"#;
        let output = vec![2.0863, 1.1865, 1544.6603, 1.2670];
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
 POSITION                   FORCE-CONSTANT FOR ION    1 DIRECTION 1 (eV/Angst/Angst)
 -----------------------------------------------------------------------------------
      0.00000      0.00000      0.00000      -330.030675  -2829.327436  -2829.327436
      0.25000      0.25000      0.25000       330.030675   2829.327436   2829.327436
 -----------------------------------------------------------------------------------
--
 POSITION                   FORCE-CONSTANT FOR ION    2 DIRECTION 1 (eV/Angst/Angst)
 -----------------------------------------------------------------------------------
      0.00000      0.00000      0.00000       329.868719   2829.280221   2829.280221
      0.25000      0.25000      0.25000      -329.868719  -2829.280221  -2829.280221
 -----------------------------------------------------------------------------------
--
 POSITION                                       TOTAL-FORCE (eV/Angst)
 -----------------------------------------------------------------------------------
      3.89220      4.01520      4.00000        -0.930834     -0.563415      0.000000
      3.00000      2.48290      4.00000        -0.006828      0.527001      0.000000
      2.12280      4.01520      4.00000         0.458533     -0.304111      0.000000
      3.00000      3.50000      4.00000         0.479129      0.340525      0.000000
 -----------------------------------------------------------------------------------
 POSITION                                       TOTAL-FORCE (eV/Angst)
 -----------------------------------------------------------------------------------
      3.86220      4.01520      4.00000         0.089245     -0.065055      0.000000
      3.00000      2.48290      4.00000         0.007618      0.545925      0.000000
      2.12280      4.01520      4.00000         0.417195     -0.352508      0.000000
      3.00000      3.50000      4.00000        -0.514057     -0.128362      0.000000
 -----------------------------------------------------------------------------------
--
"#;
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
 old parameters found on file WAVECAR:
      direct lattice vectors                 reciprocal lattice vectors
     4.001368000  0.000000000  0.000000000     0.249914529  0.000000000  0.000000000
     0.000000000  4.001368000  0.000000000     0.000000000  0.249914529  0.000000000
     0.000000000  0.000000000  4.215744000     0.000000000  0.000000000  0.237206054
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
        let output = vec![Some(vec![42.0005098f64])];
        assert_eq!(Outcar::parse_magmoms(&input), output);


        let input = r#"
 total energy-change (2. order) :-0.5897058E-05  (-0.8072299E-08)
 number of electron     309.9999998 magnetization      42.0005098 42.0005098 42.0005098
 augmentation part       88.5937960 magnetization      26.8073410 26.8073410 26.8073410
......
  free  energy   TOTEN  =      -391.79003630 eV
  energy  without entropy=     -391.77828290  energy(sigma->0) =     -391.78611850
"#;
        let output = vec![Some(vec![42.0005098f64; 3])];
        assert_eq!(Outcar::parse_magmoms(&input), output);


        let input = r#"
 total energy-change (2. order) :-0.5897058E-05  (-0.8072299E-08)
 number of electron     309.9999998 magnetization
 augmentation part       88.5937960 magnetization
......
  free  energy   TOTEN  =      -391.79003630 eV
  energy  without entropy=     -391.77828290  energy(sigma->0) =     -391.78611850
"#;
        let output = vec![None];
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
        let output = vec![Some(vec![42.0005098f64; 3]); 3];
        assert_eq!(Outcar::parse_magmoms(&input), output);
    }

    #[test]
    fn test_parse_ion_masses() {
        let input = r#"
   RPACOR =    1.200    partial core radius
   POMASS =   10.811; ZVAL   =    3.000    mass and valenz
   RCORE  =    1.700    outmost cutoff radius
--
   RPACOR =    1.200    partial core radius
   POMASS =   14.001; ZVAL   =    5.000    mass and valenz
   RCORE  =    1.500    outmost cutoff radius
--
   RPACOR =    1.200    partial core radius
   POMASS =   12.011; ZVAL   =    4.000    mass and valenz
   RCORE  =    1.500    outmost cutoff radius
--
   RPACOR =    1.800    partial core radius
   POMASS =   22.990; ZVAL   =    7.000    mass and valenz
   RCORE  =    2.200    outmost cutoff radius
--
  Mass of Ions in am
   POMASS =  10.81 14.00 12.01 22.99
  Ionic Valenz
--
   support grid    NGXF=   224 NGYF=  224 NGZF=  300
   ions per type =              18  18 108   1
 NGX,Y,Z   is equivalent  to a cutoff of  12.40, 12.40, 12.47 a.u."#;

        let output = vec![10.811; 18].into_iter()
            .chain(vec![14.001; 18].into_iter())
            .chain(vec![12.011; 108].into_iter())
            .chain(vec![22.990].into_iter())
            .collect::<Vec<_>>();

        assert_eq!(Outcar::parse_ion_masses(&input), output);
    }

    #[test]
    fn test_parse_dof() {
        let input = r#"
   Step               POTIM =    1.4999999999999999E-002
   Degrees of freedom DOF   =           3
  LATTYP: Found a simple orthorhombic cell. "#;
        let output = Some(3i32);
        assert_eq!(Outcar::_parse_dof(&input), output);
    }

    #[test]
    fn test_parse_single_vibmode() {
        let input = r#"
   2 f  =  108.545068 THz   682.008775 2PiTHz 3620.673620 cm-1   448.906478 meV
             X         Y         Z           dx          dy          dz
      3.877200  4.015200  4.000000     0.577374    0.346813    0.000001
      3.000000  2.482900  4.000000    -0.016790    0.000464    0.000000
      2.122800  4.015200  4.000000     0.577337   -0.346802   -0.000001
      3.000000  3.500000  4.000000    -0.304117   -0.000127   -0.000000 "#;
        let output = Viberation::new(3620.673620f64,
                                     vec![[ 0.577374,   0.346813,   0.000001],
                                          [-0.016790,   0.000464,   0.000000],
                                          [ 0.577337,  -0.346802,  -0.000001],
                                          [-0.304117,  -0.000127,  -0.000000]], false);

        assert_eq!(Outcar::_parse_single_vibmode(&input), output);

        let input = r#"
  10 f/i=    0.022552 THz     0.141700 2PiTHz    0.752260 cm-1     0.093268 meV
             X         Y         Z           dx          dy          dz
      3.877200  4.015200  4.000000    -0.000213    0.242665   -0.002062
      3.000000  2.482900  4.000000    -0.000118    0.242678   -0.002057
      2.122800  4.015200  4.000000    -0.000027    0.242662   -0.002062
      3.000000  3.500000  4.000000    -0.000445    0.907339   -0.007730 "#;

        let output = Viberation::new(0.752260f64,
                                     vec![[-0.000213,   0.242665,  -0.002062],
                                          [-0.000118,   0.242678,  -0.002057],
                                          [-0.000027,   0.242662,  -0.002062],
                                          [-0.000445,   0.907339,  -0.007730]], true);
        assert_eq!(Outcar::_parse_single_vibmode(&input), output);
    }

    #[test]
    fn test_parse_viberations() {
        let input = r#"
   RPACOR =    0.000    partial core radius
   POMASS =    1.000; ZVAL   =    1.000    mass and valenz
   RCORE  =    1.100    outmost cutoff radius
--
   RPACOR =    1.200    partial core radius
   POMASS =   14.001; ZVAL   =    5.000    mass and valenz
   RCORE  =    1.500    outmost cutoff radius
--
  Mass of Ions in am
   POMASS =   1.00 14.00
  Ionic Valenz

   support grid    NGXF=    60 NGYF=   72 NGZF=   80
   ions per type =               3   1
 NGX,Y,Z   is equivalent  to a cutoff of   8.31,  8.55,  8.31 a.u.

   Step               POTIM =    1.4999999999999999E-002
   Degrees of freedom DOF   =           3
  LATTYP: Found a simple orthorhombic cell.

 Eigenvectors and eigenvalues of the dynamical matrix
 ----------------------------------------------------


   1 f  =  108.762017 THz   683.371905 2PiTHz 3627.910256 cm-1   449.803706 meV
             X         Y         Z           dx          dy          dz
      3.877200  4.015200  4.000000    -0.351753   -0.188283   -0.000001
      3.000000  2.482900  4.000000    -0.000006   -0.766624    0.000001
      2.122800  4.015200  4.000000     0.352227   -0.188565   -0.000001
      3.000000  3.500000  4.000000    -0.000124    0.305756    0.000000

   2 f  =  108.545068 THz   682.008775 2PiTHz 3620.673620 cm-1   448.906478 meV
             X         Y         Z           dx          dy          dz
      3.877200  4.015200  4.000000     0.577374    0.346813    0.000001
      3.000000  2.482900  4.000000    -0.016790    0.000464    0.000000
      2.122800  4.015200  4.000000     0.577337   -0.346802   -0.000001
      3.000000  3.500000  4.000000    -0.304117   -0.000127   -0.000000

  10 f/i=    0.022552 THz     0.141700 2PiTHz    0.752260 cm-1     0.093268 meV
             X         Y         Z           dx          dy          dz
      3.877200  4.015200  4.000000    -0.000213    0.242665   -0.002062
      3.000000  2.482900  4.000000    -0.000118    0.242678   -0.002057
      2.122800  4.015200  4.000000    -0.000027    0.242662   -0.002062
      3.000000  3.500000  4.000000    -0.000445    0.907339   -0.007730

 Eigenvectors after division by SQRT(mass)

 Eigenvectors and eigenvalues of the dynamical matrix
 ----------------------------------------------------


   1 f  =  108.762017 THz   683.371905 2PiTHz 3627.910256 cm-1   449.803706 meV
             X         Y         Z           dx          dy          dz
      3.877200  4.015200  4.000000    -0.351753   -0.188283   -0.000001
      3.000000  2.482900  4.000000    -0.000006   -0.766624    0.000001
      2.122800  4.015200  4.000000     0.352227   -0.188565   -0.000001
      3.000000  3.500000  4.000000    -0.000033    0.081714    0.000000

   2 f  =  108.545068 THz   682.008775 2PiTHz 3620.673620 cm-1   448.906478 meV
             X         Y         Z           dx          dy          dz
      3.877200  4.015200  4.000000     0.577374    0.346813    0.000001
      3.000000  2.482900  4.000000    -0.016790    0.000464    0.000000
      2.122800  4.015200  4.000000     0.577337   -0.346802   -0.000001
      3.000000  3.500000  4.000000    -0.081276   -0.000034   -0.000000

  10 f/i=    0.022552 THz     0.141700 2PiTHz    0.752260 cm-1     0.093268 meV
             X         Y         Z           dx          dy          dz
      3.877200  4.015200  4.000000    -0.000213    0.242665   -0.002062
      3.000000  2.482900  4.000000    -0.000118    0.242678   -0.002057
      2.122800  4.015200  4.000000    -0.000027    0.242662   -0.002062
      3.000000  3.500000  4.000000    -0.000119    0.242488   -0.002066

 Finite differences POTIM=   1.4999999999999999E-002
  LATTYP: Found a simple orthorhombic cell.
"#;
        let masses_sqrt = vec![1.0f64, 1.0, 1.0, 14.001]
            .into_iter()
            .map(|x| x.sqrt())
            .collect::<Vec<_>>();
        let freqs = vec![3627.910256, 3620.673620, 0.752260];
        let dxdydzs =
            vec![vec![[-0.351753/masses_sqrt[0],  -0.188283/masses_sqrt[0],  -0.000001/masses_sqrt[0]],
                      [-0.000006/masses_sqrt[1],  -0.766624/masses_sqrt[1],   0.000001/masses_sqrt[1]],
                      [ 0.352227/masses_sqrt[2],  -0.188565/masses_sqrt[2],  -0.000001/masses_sqrt[2]],
                      [-0.000124/masses_sqrt[3],   0.305756/masses_sqrt[3],   0.000000/masses_sqrt[3]]],
                 vec![[ 0.577374/masses_sqrt[0],   0.346813/masses_sqrt[0],   0.000001/masses_sqrt[0]],
                      [-0.016790/masses_sqrt[1],   0.000464/masses_sqrt[1],   0.000000/masses_sqrt[1]],
                      [ 0.577337/masses_sqrt[2],  -0.346802/masses_sqrt[2],  -0.000001/masses_sqrt[2]],
                      [-0.304117/masses_sqrt[3],  -0.000127/masses_sqrt[3],  -0.000000/masses_sqrt[3]]],
                 vec![[-0.000213/masses_sqrt[0],   0.242665/masses_sqrt[0],  -0.002062/masses_sqrt[0]],
                      [-0.000118/masses_sqrt[1],   0.242678/masses_sqrt[1],  -0.002057/masses_sqrt[1]],
                      [-0.000027/masses_sqrt[2],   0.242662/masses_sqrt[2],  -0.002062/masses_sqrt[2]],
                      [-0.000445/masses_sqrt[3],   0.907339/masses_sqrt[3],  -0.007730/masses_sqrt[3]]]];
        let is_imagines = vec![false, false, true];

        let output = Some(
            freqs.into_iter()
                 .zip(dxdydzs.into_iter())
                 .zip(is_imagines.into_iter())
                 .map(|((f, d), im)| Viberation::new(f, d, im))
                 .collect::<Vec<_>>()
        );

        assert_eq!(Outcar::parse_viberations(&input), output);


        let input = r#"
   RPACOR =    0.000    partial core radius
   POMASS =    1.000; ZVAL   =    1.000    mass and valenz
   RCORE  =    1.100    outmost cutoff radius
--
   RPACOR =    1.200    partial core radius
   POMASS =   14.001; ZVAL   =    5.000    mass and valenz
   RCORE  =    1.500    outmost cutoff radius
--
  Mass of Ions in am
   POMASS =   1.00 14.00
  Ionic Valenz

   support grid    NGXF=    60 NGYF=   72 NGZF=   80
   ions per type =               3   1
 NGX,Y,Z   is equivalent  to a cutoff of   8.31,  8.55,  8.31 a.u.

   Step               POTIM =    1.4999999999999999E-002
   Degrees of freedom DOF   =           3
  LATTYP: Found a simple orthorhombic cell.
"#;
        let output = None;
        assert_eq!(Outcar::parse_viberations(&input), output);
    }
}
