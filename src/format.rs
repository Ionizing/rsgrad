use std::fmt;
use std::path::{
    Path,
    PathBuf,
};
use std::fs;
use std::io;
use std::io::Write;

use itertools::multizip;
use colored::Colorize;
use vasp_poscar::{self, Poscar};
use crate::outcar::{
    Outcar,
    IonicIteration,
    Mat33,
    MatX3,
};

#[derive(Clone)]
pub struct IonicIterationsFormat {
    _data            : Vec<IonicIteration>,

    print_energy     : bool,
    print_energyz    : bool,
    print_log10de    : bool,
    print_favg       : bool,
    print_fmax       : bool,
    print_fmax_axis  : bool,
    print_fmax_index : bool,
    print_nscf       : bool,
    print_time_usage : bool,
    print_magmom     : bool,
    print_volume     : bool,
}

impl From<Vec<IonicIteration>> for IonicIterationsFormat {
    fn from(data: Vec<IonicIteration>) -> Self {
        Self {
            _data            : data,
            print_energy     : false,
            print_energyz    : true,
            print_log10de    : false,
            print_favg       : true,
            print_fmax       : true,
            print_fmax_axis  : false,
            print_fmax_index : false,
            print_nscf       : true,
            print_time_usage : true,
            print_magmom     : true,
            print_volume     : false,
        }
    }
}

macro_rules! impl_builder_item {
    ($t: tt) => {
        pub fn $t(mut self, arg: bool) -> Self {
            self.$t = arg;
            self
        }
    };
}

// Use non-consuming builder pattern
impl IonicIterationsFormat {
    impl_builder_item!(print_energy);
    impl_builder_item!(print_energyz);
    impl_builder_item!(print_log10de);
    impl_builder_item!(print_favg);
    impl_builder_item!(print_fmax);
    impl_builder_item!(print_fmax_axis);
    impl_builder_item!(print_fmax_index);
    impl_builder_item!(print_nscf);
    impl_builder_item!(print_time_usage);
    impl_builder_item!(print_magmom);
    impl_builder_item!(print_volume);
}

impl fmt::Display for IonicIterationsFormat {
    fn fmt(&self, f: &mut fmt::Formatter<'_>) -> fmt::Result {
        let mut ce: f64 = 0.0;

        // Prepare Header
        let mut header = "  #Step".to_owned();
        header += if self.print_energy     { "    TOTEN/eV" } else { "" };
        header += if self.print_energyz    { "  TOTEN_z/eV" } else { "" };
        header += if self.print_log10de    { " LgdE" }        else { "" };
        header += if self.print_favg       { "   Favg" }      else { "" };
        header += if self.print_fmax       { "   Fmax" }      else { "" };
        header += if self.print_fmax_index { " idx" }         else { "" };
        header += if self.print_fmax_axis  { "  " }           else { "" };
        header += if self.print_nscf       { " #SCF" }        else { "" };
        header += if self.print_time_usage { " Time/m" }      else { "" };
        header += if self.print_volume     { "   Vol/A3" }    else { "" };
        header += if self.print_magmom     { " Mag/muB" }     else { "" };
        writeln!(f, "{}", header.bright_green())?;

        for (i, it) in self._data.iter().enumerate() {
            let mut line = format!("{:7}", i+1);

            let de = self._data[i].toten_z - ce;
            ce = self._data[i].toten_z;
            if self.print_energy  { line += &format!(" {:11.5}", it.toten); }
            if self.print_energyz { line += &format!(" {:11.5}", it.toten_z).bright_green().to_string(); }
            if self.print_log10de { line += &format!(" {:4.1}", de.abs().log10()); }

            let fsize = it.forces.iter()
                                   .map(|f| (f[0]*f[0] + f[1]*f[1] * f[2]*f[2]).sqrt())
                                   .collect::<Vec<_>>();

            if self.print_favg {
                line += &format!(" {:6.3}", fsize.iter().sum::<f64>() / it.forces.len() as f64);
            }

            let (fmax_ind, fmax) = fsize.into_iter()
                                        .enumerate()
                                        .fold((0, 0.0), |mut acc, (i, f)|{
                                            if acc.1 < f {
                                                acc.1 = f;
                                                acc.0 = i;
                                            }
                                            acc
                                        });
            let fmaxis = match it.forces[fmax_ind]
                .iter()
                .enumerate()
                .fold((0, 0.0), |mut acc, (i, f)| {
                    if acc.1 < f.abs() {
                        acc.1 = f.abs();
                        acc.0 = i;
                    }
                    acc
                }) {
                    (0, _) => "X",
                    (1, _) => "Y",
                    (2, _) => "Z",
                    _ => unreachable!("Invalid Fmax Axis here")
                };

            if self.print_fmax       { line += &format!(" {:6.3}", fmax).bright_green().to_string(); }
            if self.print_fmax_index { line += &format!(" {:3}", fmax_ind+1); }
            if self.print_fmax_axis  { line += &format!(" {:1}", fmaxis); }
            if self.print_nscf       { line += &format!(" {:4}", it.nscf).bright_yellow().to_string(); }
            if self.print_time_usage { line += &format!(" {:6.2}", it.cputime/60.0); }

            if self.print_volume {
                let volume = {
                    let c = it.cell;

                    // |00 01 02|
                    // |10 11 12|
                    // |20 21 22|

                    c[0][0] * (c[1][1] * c[2][2] - c[2][1] * c[1][2])
                        - c[0][1] * (c[1][0] * c[2][2] - c[1][2] * c[2][0])
                        + c[0][2] * (c[1][0] * c[2][1] - c[1][1] * c[2][0])
                };
                line += &format!(" {:8.1}", volume);
            }

            if self.print_magmom {
                if let Some(mag) = &it.magmom {
                    line += &mag.iter()
                                .map(|n| format!(" {:7.3}", n))
                                .collect::<Vec<_>>()
                                .join("");
                } else { line += "   NoMag"; }
            }

            writeln!(f, "{}", line)?;
        }
        Ok(())
    }

}


#[derive(Clone)]
pub struct Structure {
    pub cell          : Mat33<f64>,
    pub ion_types     : Vec<String>,
    pub ions_per_type : Vec<i32>,
    pub car_pos       : MatX3<f64>,
    pub frac_pos      : MatX3<f64>,
}


impl Outcar {
    pub fn get_structure_cloned(&self, i: usize) -> Structure {
        let cell     = self.ion_iters[i].cell.clone();
        let car_pos  = self.ion_iters[i].positions.clone();
        let frac_pos = _car_to_frac(&cell, &car_pos);

        Structure {
            cell,
            ion_types: self.ion_types.clone(),
            ions_per_type: self.ions_per_type.clone(),
            car_pos,
            frac_pos,
        }
    }

}

fn _car_to_frac(cell: &Mat33<f64>, carpos: &MatX3<f64>) -> MatX3<f64> {
    let convmat = _calc_inv_3x3(cell);
    carpos.iter()
          .map(|v| {
              [
                  //                    [B11, B12, B13]
                  // carpos [x, y, z] x [B21, B22, B23]
                  //                    [B31, B32, B33]
                  convmat[0][0] * v[0] + convmat[1][0] * v[1] + convmat[2][0] * v[2],
                  convmat[0][1] * v[0] + convmat[1][1] * v[1] + convmat[2][1] * v[2],
                  convmat[0][2] * v[0] + convmat[1][2] * v[1] + convmat[2][2] * v[2],
              ]
          }).collect()
}

fn _calc_inv_3x3(cell: &Mat33<f64>) -> Mat33<f64> {
    let a = cell[0][0];
    let b = cell[0][1];
    let c = cell[0][2];
    let d = cell[1][0];
    let e = cell[1][1];
    let f = cell[1][2];
    let g = cell[2][0];
    let h = cell[2][1];
    let i = cell[2][2];

    // [a b c]
    // [d e f]
    // [g h i]

    let det = a*(e*i - f*h) - b*(d*i - f*g) + c*(d*h - e*g);
    assert!(det.abs() > 1.0e-5, "Cell volume too small, check your lattice.");
    let detdiv = 1.0 / det;

    [[(e*i - f*h) * detdiv, (c*h - b*i) * detdiv, (b*f - c*e) * detdiv],
     [(f*g - d*i) * detdiv, (a*i - c*g) * detdiv, (c*d - a*f) * detdiv],
     [(d*h - e*g) * detdiv, (b*g - a*h) * detdiv, (a*e - b*d) * detdiv]]
}

#[derive(Clone)]
pub struct Trajectory(Vec<Structure>);

impl Trajectory {
    pub fn save_as_xdatcar(self, path: &(impl AsRef<Path> + ?Sized)) -> io::Result<()> {
        let mut f = fs::OpenOptions::new()
            .create(true)
            .truncate(true)
            .write(true)
            .open(path)?;

        for (i, v) in self.0.into_iter().enumerate() {
            //
            // ------
            // Generated by rsgrad
            //    1.000000
            //    [ax, ay, az]
            //    [bx, by, bz]
            //    [cx, cy, cz]
            //    H
            //    1
            // Direct configuration=     1
            //  0.00000000 0.00000000 0.00000000
            // Generated by rsgrad
            //    1.000000
            //    [ax, ay, az]
            //    [bx, by, bz]
            //    [cx, cy, cz]
            //    H
            //    1
            // Direct configuration=     1
            //  0.00000000 0.00000000 0.00000000
            // ...
            // ...
            // ------
            writeln!(f, "Generated by rsgrad")?;
            writeln!(f, "{:15.9}", 1.0)?;
            for row in v.cell.iter() {
                writeln!(f, " {:12.6}{:12.6}{:12.6}", row[0], row[1], row[2])?;
            }

            for elem in v.ion_types {
                write!(f, "{:>4}", elem)?;
            }
            writeln!(f, "")?;
            for nelm in v.ions_per_type {
                write!(f, "{:>4}", nelm)?;
            }
            writeln!(f, "")?;

            writeln!(f, "Direct configuration={:6}", i+1)?;
            for row in v.frac_pos {
                writeln!(f, " {:15.9} {:15.9} {:15.9}", row[0], row[1], row[2])?;
            }
        }
        Ok(())
    }

    pub fn save_as_seperated_poscars(self, path: &(impl AsRef<Path> + ?Sized)) -> io::Result<()> {
        let mut fname = PathBuf::new();
        fname.push(path);
        if !fname.is_dir() {
            fs::create_dir_all(&fname)?;
        }

        for (i, v) in self.0.into_iter().enumerate() {
            fname.push(&format!("POSCAR_{:05}.vasp", i+1));
            let mut f = fs::OpenOptions::new()
                .create(true)
                .truncate(true)
                .write(true)
                .open(&fname)?;
            write!(f, "{:15.9}", Poscar::from(v))?;
            fname.pop();
        }
        Ok(())
    }

    pub fn _save_into_seperated_dirs(&self, path: &(impl AsRef<Path> + ?Sized)) -> io::Result<()> {
        todo!();
    }
}


impl From<Outcar> for Trajectory {
    fn from(o: Outcar) -> Self {
        let len = o.ion_iters.len();
        // todo!();
        Self (
            multizip((o.ion_iters, vec![o.ion_types; len], vec![o.ions_per_type; len]))
                .map(|(ii, it, ipt)| -> Structure {
                    let cell = ii.cell;
                    let car_pos = ii.positions;
                    let frac_pos = _car_to_frac(&cell, &car_pos);
                    Structure {
                        cell,
                        ion_types: it,
                        ions_per_type: ipt,
                        car_pos,
                        frac_pos,
                    }
                })
                .collect()
        )
    }
}


impl From<Structure> for Poscar {
    fn from(s: Structure) -> Self {
        vasp_poscar::Builder::new()
            .comment("Generated by rsgrad")
            .scale(vasp_poscar::ScaleLine::Factor(1.0))
            .lattice_vectors(&s.cell)
            .group_symbols(s.ion_types)
            .group_counts(s.ions_per_type.into_iter()
                          .map(|x| x as usize)
                          .collect::<Vec<usize>>())
            .positions(
                vasp_poscar::Coords::Cart(s.car_pos)
            )
            .build()
            .unwrap()
    }
}


impl Structure {
    pub fn save_as_poscar(self, path: &(impl AsRef<Path> + ?Sized)) -> io::Result<()> {
        let mut f = fs::OpenOptions::new()
            .create(true)
            .truncate(true)
            .write(true)
            .open(path)?;
        write!(f, "{:.9}", Poscar::from(self))
    }
}


#[cfg(test)]
mod tests {
    use super::*;

    #[test]
    fn test_structure_to_poscar() {
        let s = Structure {
            cell: [[5.0, 0.0, 0.0], [0.0, 5.0, 0.0], [0.0, 0.0, 5.0]],
            ion_types: vec!["H".to_string()],
            ions_per_type: vec![1],
            car_pos: vec![[0.0, 0.0, 0.0]],
            frac_pos: vec![[0.0, 0.0, 0.0]],
        };
        println!("{:15.9}", Poscar::from(s.clone()));
        assert_eq!(r#"Generated by rsgrad
      1.000000000
        5.000000000     0.000000000     0.000000000
        0.000000000     5.000000000     0.000000000
        0.000000000     0.000000000     5.000000000
   H
   1
Cartesian
      0.000000000     0.000000000     0.000000000
"#, format!("{:15.9}", Poscar::from(s)));
    }

    #[test]
    fn test_calc_inv_3x3() {
        let cell = [[1.0, 2.0, 3.0],
                    [0.0, 1.0, 4.0],
                    [5.0, 6.0, 0.0]];
        assert_eq!(_calc_inv_3x3(&cell), [[-24.0,  18.0,  5.0],
                                          [ 20.0, -15.0, -4.0],
                                          [ -5.0,   4.0,  1.0]]);
    }

    #[test]
    #[should_panic]
    fn test_inv_3x3_singular() {
        let cell = [[1.0, 2.0, 3.0],
                    [4.0, 5.0, 6.0],
                    [7.0, 8.0, 9.0]];
        let _ = _calc_inv_3x3(&cell);
    }

    #[test]
    fn test_car_to_frac() {
        let cell = [[1.0, 2.0, 3.0],
                    [0.0, 1.0, 4.0],
                    [5.0, 6.0, 0.0]];
        let car = vec![[0.0, 0.0, 0.0],
                       [0.5, 1.0, 1.5],
                       [0.0, 0.5, 2.0],
                       [2.5, 3.0, 0.0],
                       [11.0, 45.0, 14.0]];
        let frac = vec![[0.0, 0.0, 0.0],
                        [0.5, 0.0, 0.0],
                        [0.0, 0.5, 0.0],
                        [0.0, 0.0, 0.5],
                        [566.0, -421.0, -111.0]];
        assert_eq!(frac, _car_to_frac(&cell, &car));
    }
}
