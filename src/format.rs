use std::fmt;

use crate::outcar::IonicIteration;

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
        pub fn $t(&mut self, arg: bool) -> &mut Self {
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
        let mut de: f64 = 0.0;

        for (i, it) in self._data.iter().enumerate() {
            let mut line = format!("{:7}", i+1);

            de -= self._data[i].toten_z;
            if self.print_energy  { line += &format!(" {:11.5}", it.toten); }
            if self.print_energyz { line += &format!(" {:11.5}", it.toten_z); }
            if self.print_log10de { line += &format!(" {:4.1}", de.log10()); }

            let fsize = it.forces.iter()
                                   .map(|f| (f[0]*f[0] + f[1]*f[1] * f[2]*f[2]).sqrt())
                                   .collect::<Vec<_>>();

            if self.print_favg {
                line += &format!(" {:6.3}", fsize.iter().sum::<f64>() / it.forces.len() as f64);
            }

            let (fmax_ind, fmax) = fsize.into_iter()
                                        .enumerate()
                                        .fold((0, 0.0), |mut acc, (i, f)|{
                                            if acc.1 > f {
                                                acc.1 = f;
                                                acc.0 = i;
                                            }
                                            acc
                                        });

            if self.print_fmax       { line += &format!(" {:6.3}", fmax); }
            if self.print_fmax_index { line += &format!(" {:3}", fmax_ind); }
            if self.print_nscf       { line += &format!(" {:3}", it.nscf); }
            if self.print_time_usage { line += &format!(" {:6.2}", it.cputime/60.0); }

            if self.print_volume {
                let volume = {
                    let c = it.cell;

                    // |00 01 02|
                    // |10 11 12|
                    // |20 21 22|

                    c[0][0] * (c[1][1] * c[2][2] - c[2][1] * c[1][2]) +
                        c[0][1] * (c[1][0] * c[2][2] - c[1][2] * c[2][0]) +
                        c[0][2] * (c[1][0] * c[2][1] - c[1][1] * c[2][0])
                };
                line += &format!(" {}", volume);
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




// #[cfg(test)]
// mod tests {
//     use super::*;
//     #[test]
//     fn test_into() {
//     }
// }
