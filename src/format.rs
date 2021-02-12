use std::fmt;

use itertools::multizip;
use colored::Colorize;
use vasp_poscar::{
    self, Poscar
};
use crate::outcar::{
    Outcar,
    IonicIteration,
    Mat33,
    MatX3,
};

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


pub struct Structure {
    pub cell          : Mat33<f64>,
    pub ion_types     : Vec<String>,
    pub ions_per_type : Vec<i32>,
    pub positions     : MatX3<f64>,
}


impl Outcar {
    pub fn get_structure_cloned(&self, i: usize) -> Structure {
        Structure {
            cell: self.ion_iters[i].cell.clone(),
            ion_types: self.ion_types.clone(),
            ions_per_type: self.ions_per_type.clone(),
            positions: self.ion_iters[i].positions.clone(),
        }
    }
}


impl From<Outcar> for Vec<Structure> {
    fn from(o: Outcar) -> Self {
        let len = o.ion_iters.len();
        multizip((o.ion_iters, vec![o.ion_types; len], vec![o.ions_per_type; len]))
            .map(|(ii, it, ipt)| -> Structure {
                Structure {
                    cell: ii.cell,
                    ion_types: it,
                    ions_per_type: ipt,
                    positions: ii.positions,
                }
            })
            .collect()
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
                vasp_poscar::Coords::Cart(s.positions)
            )
            .build()
            .unwrap()
    }
}
