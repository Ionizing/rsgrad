use std::{
    fs,
    path::Path,
};
use regex::Regex;
use ndarray::{
    s,
    Array5,
    Array6,
};
use anyhow::{
    Result,
    Context,
};
use crate::types::{
    Vector,
    Matrix,
    Cube,
};
use crate::vasp_parsers::kpoints::Kpoints;

#[derive(Clone)]
pub struct ProjectedDOS {
    pub nions:       u32,
    pub nspin:       u32,
    pub nkpoints:    u32,
    pub nbands:      u32,
    pub lsorbit:     bool,
    pub nlm:         Vec<String>,
    pub eigvals:     Cube<f64>,             // [ispin, ikpoint, iband]
    pub occupations: Cube<f64>,
    pub projected:   Array5<f64>,     // [[ispin, ikpoint, iband], [iion, iorbit]]
}

#[derive(Clone)]
pub struct Procar {
    pub kpoints: Kpoints,
    pub pdos:    ProjectedDOS,
}

impl Procar {
    pub fn from_file(path: impl AsRef<Path>) -> Result<Self> {
        let txt = fs::read_to_string(path)?;

        let nspin                       = Self::parse_nspin(&txt)?;
        let lsorbit                     = Self::parse_lsorbit(&txt)?;
        let (nkpoints, nbands, nions)   = Self::parse_kbi(&txt)?;
        let nlm                         = Self::parse_nlm(&txt)?;
        let projected                   = Self::parse_projections(&txt, nspin, nkpoints, nbands, nions, nlm.len(), lsorbit)?;
        let (eigvals, occupations)      = Self::parse_bandlevel(&txt, nspin, nkpoints, nbands)?;

        let (kpointlist, weights) = Self::parse_kpoints(&txt, nkpoints)?;

        let kpoints = Kpoints {
            nkpoints: nkpoints as u32,
            kpointlist,
            weights,
        };

        let pdos = ProjectedDOS {
            nions: nions as u32,
            nspin: nspin as u32,
            nkpoints: nkpoints as u32,
            nbands: nbands as u32,
            lsorbit,
            nlm,
            eigvals,
            occupations,
            projected,
        };

        Ok(Self{
            kpoints,
            pdos,
        })
    }

    fn parse_nspin(txt: &str) -> Result<usize> {
        Ok(
            Regex::new(r"(?m)^ k-point \s+1 :")?
            .find_iter(txt)
            .count()
          )
    }

    fn parse_lsorbit(txt: &str) -> Result<bool> {
        let range = Regex::new(r"(?m)^band ")?
            .find_iter(txt)
            .take(2)
            .map(|m| m.start())
            .collect::<Vec<_>>();
        let range = range[0] .. range[1];

        Ok(
            Regex::new(r"(?m)^tot ")?
            .find_iter(&txt[range])
            .count() == 4
          )
    }

    fn parse_kbi(txt: &str) -> Result<(usize, usize, usize)> {
        let caps = Regex::new(r"(?m)^# of k-points:\s*(\d+)\s+# of bands:\s*(\d+)\s+# of ions:\s*(\d+)")?
            .captures(txt)
            .context("Cannot find metadata of PROCAR.")?;
        let nkpoints = caps.get(1).unwrap().as_str().parse()?;
        let nbands   = caps.get(2).unwrap().as_str().parse()?;
        let nions    = caps.get(3).unwrap().as_str().parse()?;
        Ok((nkpoints, nbands, nions))
    }

    fn parse_nlm(txt: &str) -> Result<Vec<String>> {
        let mut ret = txt.lines()
            .nth(7).context("Cannot find orbit tags from PROCAR.")?
            .split_whitespace()
            .skip(1)
            .map(str::to_string)
            .collect::<Vec<_>>();
        ret.pop();
        Ok(ret)
    }

    fn parse_projections(txt: &str, nspin: usize, nkpoints: usize,
                                    nbands: usize, nions: usize, nnlm: usize, lsorbit: bool) -> Result<Array5<f64>> {
        let mut pdos = Vec::<Matrix<f64>>::with_capacity(16);
        for p in Regex::new(r"(?m)^band")?.find_iter(txt) {
            pdos.push(Procar::parse_projections_aux(&txt[p.start()..], nions, lsorbit)?);
        }

        let pdos = pdos.into_iter().flatten().collect::<Vec<_>>();
        let pdos = if !lsorbit {
            Array5::from_shape_vec((nspin, nkpoints, nbands, nions, nnlm), pdos)?
        } else {
            let mut pdos6 = Array6::from_shape_vec((1, nkpoints, nbands, 4, nions, nnlm), pdos)?;
            pdos6.swap_axes(0, 3);
            let pdos = pdos6.into_iter().collect::<Vec<_>>();
            Array5::from_shape_vec((4, nkpoints, nbands, nions, nnlm), pdos)?
        };

        Ok(pdos)
    }

    fn parse_projections_aux(txt: &str, nions: usize, lsorbit: bool) -> Result<Matrix<f64>> {
        let reg = Regex::new(r"^\s+\d")?;
        let nspinor = if lsorbit { 4usize } else { 1 };

        let pdos = txt.lines()
            .filter(|l| reg.is_match(l))
            .take(nspinor * nions)
            .map(|l| {
                l.split_whitespace()
                 .map(|t| t.parse::<f64>().unwrap())
            })
            .flatten()
            .collect::<Vec<_>>();

        let ncol = pdos.len() / nions / nspinor;
        let pdos = Matrix::from_shape_vec((nspinor * nions, ncol), pdos)?;
        let pdos = pdos.slice(s!(.., 1..-1)).into_owned();

        Ok(pdos)
    }

    fn parse_kpoints(txt: &str, nkpoints: usize) -> Result<(Matrix<f64>, Vector<f64>)> {
        let mut kpoints = Vec::<f64>::with_capacity(16);
        let mut weights = Vec::<f64>::with_capacity(16);

        Regex::new(r"(?m)^ k-point .{5,6}:   ([- ]\d\.\d{8})([- ]\d\.\d{8})([- ]\d\.\d{8})     weight = (\d\.\d{8})$")?
            .captures_iter(txt)
            .take(nkpoints)
            .for_each(|m| {
                kpoints.extend(&[
                    m.get(1).unwrap().as_str().trim().parse::<f64>().unwrap(),
                    m.get(2).unwrap().as_str().trim().parse::<f64>().unwrap(),
                    m.get(3).unwrap().as_str().trim().parse::<f64>().unwrap(),
                ]);
                weights.push(m.get(4).unwrap().as_str().trim().parse::<f64>().unwrap());
            });

        let kpoints = Matrix::from_shape_vec((nkpoints, 3), kpoints)?;
        let weights = weights.into();
        Ok((kpoints, weights))
    }

    fn parse_bandlevel(txt: &str, nspin: usize, nkpoints: usize, nbands: usize) -> Result<(Cube<f64>, Cube<f64>)> {
        let mut levels = Vec::<f64>::with_capacity(64);
        let mut occups = Vec::<f64>::with_capacity(64);

        Regex::new(r"(?m)^band.{4,6} # energy\s*([+-]?([0-9]+([.][0-9]*)?|[.][0-9]+)) # occ\.\s*([+-]?([0-9]+([.][0-9]*)?|[.][0-9]+))")?
            .captures_iter(&txt)
            .for_each(|m| {
                levels.push(m.get(1).unwrap().as_str().parse::<f64>().unwrap());
                occups.push(m.get(5).unwrap().as_str().parse::<f64>().unwrap());
            });

        let levels = Cube::from_shape_vec((nspin, nkpoints, nbands), levels)?;
        let occups = Cube::from_shape_vec((nspin, nkpoints, nbands), occups)?;
        Ok((levels, occups))
    }
}


#[cfg(test)]
mod test {
    use super::*;
    use ndarray::{s, arr1, arr2};

    const TXT_STD: &str = include_str!("../../tests/PROCAR");
    const TXT_SMP: &str = include_str!("../../tests/PROCAR.simple");
    const TXT_544: &str = include_str!("../../tests/PROCAR.new_format_5.4.4");
    const TXT_PHS: &str = include_str!("../../tests/PROCAR.phase");
    const TXT_NCL: &str = include_str!("../../tests/PROCAR.ncl");

    #[test]
    fn test_parse_lsorbit() -> Result<()> {
        assert_eq!(Procar::parse_lsorbit(TXT_STD)?, false);
        assert_eq!(Procar::parse_lsorbit(TXT_SMP)?, false);
        assert_eq!(Procar::parse_lsorbit(TXT_PHS)?, false);
        assert_eq!(Procar::parse_lsorbit(TXT_544)?, false);
        assert_eq!(Procar::parse_lsorbit(TXT_NCL)?, true);
        Ok(())
    }

    #[test]
    fn test_parse_nspin() -> Result<()> {
        assert_eq!(Procar::parse_nspin(TXT_STD)?, 2);
        assert_eq!(Procar::parse_nspin(TXT_SMP)?, 2);
        assert_eq!(Procar::parse_nspin(TXT_544)?, 2);
        assert_eq!(Procar::parse_nspin(TXT_PHS)?, 2);
        assert_eq!(Procar::parse_nspin(TXT_NCL)?, 1);
        Ok(())
    }

    #[test]
    fn test_parse_kbi() -> Result<()> {
        assert_eq!(Procar::parse_kbi(TXT_STD)?, (18, 49, 8));
        assert_eq!(Procar::parse_kbi(TXT_SMP)?, (10, 10, 3));
        assert_eq!(Procar::parse_kbi(TXT_544)?, (13, 20, 4));
        assert_eq!(Procar::parse_kbi(TXT_PHS)?, (60, 12, 3));
        assert_eq!(Procar::parse_kbi(TXT_NCL)?, (81, 28, 3));
        Ok(())
    }

    #[test]
    fn test_parse_nlm() -> Result<()> {
        assert_eq!(Procar::parse_nlm(TXT_STD)?,
                   vec!["s", "py", "pz", "px", "dxy", "dyz", "dz2", "dxz", "dx2"]);
        assert_eq!(Procar::parse_nlm(TXT_SMP)?,
                   vec!["s", "p", "d"]);
        assert_eq!(Procar::parse_nlm(TXT_544)?,
                   vec!["s", "py", "pz", "px", "dxy", "dyz", "dz2", "dxz", "x2-y2"]);
        assert_eq!(Procar::parse_nlm(TXT_PHS)?,
                   vec!["s", "py", "pz", "px", "dxy", "dyz", "dz2", "dxz", "dx2"]);
        assert_eq!(Procar::parse_nlm(TXT_NCL)?,
                   vec!["s", "py", "pz", "px", "dxy", "dyz", "dz2", "dxz", "x2-y2"]);
        Ok(())
    }

    #[test]
    fn test_parse_projections() -> Result<()> {
        Ok(())
    }

    #[test]
    fn test_parse_kpoints() -> Result<()> {
        let kpoints = Procar::parse_kpoints(TXT_STD, 18)?;
        assert_eq!(kpoints.0.slice(s!(0..3, 0..3)),
                   arr2(&[[0.08333333, 0.08333333, 0.08333333],
                          [0.25000000, 0.08333333, 0.08333333],
                          [0.41666667, 0.08333333, 0.08333333]]));
        assert_eq!(kpoints.0.slice(s!(-1.., 0..3)),
                   arr2(&[[0.41666667, 0.41666667, 0.41666667]]));
        assert_eq!(kpoints.1, arr1(&[0.03703704, 0.07407407, 0.07407407, 0.03703704,
                                     0.07407407, 0.07407407, 0.03703704, 0.07407407,
                                     0.07407407, 0.03703704, 0.07407407, 0.03703704,
                                     0.07407407, 0.03703704, 0.07407407, 0.03703704,
                                     0.03703704, 0.03703704]));

        let kpoints = Procar::parse_kpoints(TXT_544, 13)?;
        assert_eq!(kpoints.0.slice(s!(0..3, 0..3)),
                   arr2(&[[0.00000000, 0.00000000, 0.00000000],
                          [0.25000000, 0.00000000, 0.00000000],
                          [0.50000000, 0.00000000, 0.00000000]]));
        assert_eq!(kpoints.0.slice(s!(-1.., 0..3)),
                   arr2(&[[0.50000000, 0.50000000, 0.50000000]]));
        assert_eq!(kpoints.1, arr1(&[0.01562500, 0.09375000, 0.04687500, 0.09375000,
                                     0.18750000, 0.09375000, 0.04687500, 0.03125000,
                                     0.09375000, 0.09375000, 0.09375000, 0.09375000,
                                     0.01562500]));

        let kpoints = Procar::parse_kpoints(TXT_PHS, 60)?;
        assert_eq!(kpoints.0.slice(s!(0..3, 0..3)),
                   arr2(&[[0.06250000, 0.06250000, 0.06250000],
                          [0.18750000, 0.06250000, 0.06250000],
                          [0.31250000, 0.06250000, 0.06250000]]));
        assert_eq!(kpoints.0.slice(s!(-1.., 0..3)),
                   arr2(&[[-0.43750000, 0.43750000, 0.43750000]]));
        assert_eq!(kpoints.1.slice(s!(..3)), arr1(&[0.00390625, 0.01171875, 0.01171875]));

        let kpoints = Procar::parse_kpoints(TXT_SMP, 10)?;
        assert_eq!(kpoints.0.slice(s!(0..3, 0..3)),
                   arr2(&[[0.12500000, 0.12500000, 0.12500000],
                          [0.37500000, 0.12500000, 0.12500000],
                          [-0.3750000, 0.12500000, 0.12500000]]));
        assert_eq!(kpoints.0.slice(s!(-1.., 0..3)),
                   arr2(&[[-0.3750000, 0.37500000, 0.37500000]]));
        assert_eq!(kpoints.1.slice(s!(..3)), arr1(&[0.0312500, 0.09375000, 0.09375000]));

        let kpoints = Procar::parse_kpoints(TXT_NCL, 81)?;
        assert_eq!(kpoints.0.slice(s!(0..3, 0..3)),
                   arr2(&[[0.00000000, 0.00000000, 0.00000000],
                          [0.11111111, 0.00000000, 0.00000000],
                          [0.22222222, 0.00000000, 0.00000000]]));
        assert_eq!(kpoints.0.slice(s!(-1.., 0..3)),
                   arr2(&[[-0.11111111, -0.11111111, 0.00000000]]));
        assert_eq!(kpoints.1.slice(s!(..3)), arr1(&[0.01234568, 0.01234568, 0.01234568]));
        Ok(())
    }

    #[test]
    fn test_parse_bandlevel() -> Result<()> {
        let bandlevels = Procar::parse_bandlevel(TXT_STD, 2, 18, 49)?;
        assert_eq!(bandlevels.0[(0, 0,  0)], -51.69067761);
        assert_eq!(bandlevels.0[(0, 0,  1)], -51.68875702);
        assert_eq!(bandlevels.0[(1, 17, 48)],   9.81589913);

        assert_eq!(bandlevels.1[(0, 0, 39)], 1.0);
        assert_eq!(bandlevels.1[(0, 0, 40)], 0.0);
        assert_eq!(bandlevels.1[(1, 17, 39)], 1.0);
        assert_eq!(bandlevels.1[(1, 17, 40)], 0.0);


        let bandlevels = Procar::parse_bandlevel(TXT_SMP, 2, 10, 10)?;
        assert_eq!(bandlevels.0[(0, 0, 0)], -14.51487969);
        assert_eq!(bandlevels.0[(0, 0, 1)],   0.41588795);
        assert_eq!(bandlevels.0[(1, 9, 9)],  18.63600432);

        assert_eq!(bandlevels.1[(0, 0, 3)], 1.0);
        assert_eq!(bandlevels.1[(0, 0, 4)], 0.0);
        assert_eq!(bandlevels.1[(1, 9, 3)], 1.0);
        assert_eq!(bandlevels.1[(1, 9, 4)], 0.0);


        let bandlevels = Procar::parse_bandlevel(TXT_544, 2, 13, 20)?;
        assert_eq!(bandlevels.0[(0,  0,  0)], -12.28198856);
        assert_eq!(bandlevels.0[(0,  0,  1)], -10.87401691);
        assert_eq!(bandlevels.0[(1, 12, 19)],  12.47606922);

        assert_eq!(bandlevels.1[(0, 0,  14)], 1.0);
        assert_eq!(bandlevels.1[(0, 0,  15)], 0.99999999);
        assert_eq!(bandlevels.1[(1, 12, 14)], 0.99995168);
        assert_eq!(bandlevels.1[(1, 12, 15)], 0.99995161);

        let bandlevels = Procar::parse_bandlevel(TXT_PHS, 2, 60, 12)?;
        assert_eq!(bandlevels.0[(0,  0,  0)], -43.00085757);
        assert_eq!(bandlevels.0[(0,  0,  1)], -42.67172300);
        assert_eq!(bandlevels.0[(1, 59, 11)],  17.86374001);

        assert_eq!(bandlevels.1[(0,  0, 5)], 1.0);
        assert_eq!(bandlevels.1[(0,  0, 6)], 0.0);
        assert_eq!(bandlevels.1[(1, 59, 5)], 1.0);
        assert_eq!(bandlevels.1[(1, 59, 6)], 0.0);

        let bandlevels = Procar::parse_bandlevel(TXT_NCL, 1, 81, 28)?;
        assert_eq!(bandlevels.0[(0,  0,  0)], -16.82688427);
        assert_eq!(bandlevels.0[(0,  0,  1)], -16.81850590);
        assert_eq!(bandlevels.0[(0, 80, 27)],   2.83754252);

        assert_eq!(bandlevels.1[(0,  0, 15)], 0.00001319);
        assert_eq!(bandlevels.1[(0,  0, 16)], 0.00001256);
        assert_eq!(bandlevels.1[(0, 80, 16)], 1.00000000);
        assert_eq!(bandlevels.1[(0, 80, 17)], 0.00000002);

        Ok(())
    }

    #[test]
    fn test_parse_projected_aux() -> Result<()> {
        let txt: &str = r#"

ion      s     py     pz     px    dxy    dyz    dz2    dxz  x2-y2    tot
    1  0.060  0.000  0.000  0.000  0.000  0.000  0.000  0.000  0.000  0.060
    2  0.059  0.000  0.000  0.000  0.000  0.000  0.000  0.000  0.000  0.059
    3  0.339  0.000  0.000  0.000  0.000  0.000  0.000  0.000  0.000  0.339
    4  0.339  0.000  0.000  0.000  0.000  0.000  0.000  0.000  0.000  0.339
tot    0.796  0.000  0.000  0.000  0.000  0.000  0.000  0.000  0.000  0.796
ion          s             py             pz             px             dxy            dyz            dz2            dxz          dx2-y2
    1 -0.130  0.199  -0.000 -0.000   0.000  0.000   0.000 -0.000   0.000 -0.000   0.000 -0.000   0.000  0.000   0.000 -0.000   0.000  0.000   0.057
    2 -0.129  0.197  -0.000 -0.000  -0.000  0.000  -0.000 -0.000  -0.000  0.000  -0.000  0.000  -0.000  0.000  -0.000  0.000   0.000 -0.000   0.055
    3 -0.314  0.480  -0.001  0.001  -0.001  0.001  -0.001  0.001   0.000  0.000   0.000  0.000   0.000  0.000   0.000  0.000   0.000  0.000   0.329
    4 -0.314  0.480   0.001 -0.001   0.001 -0.001   0.001 -0.001   0.000  0.000   0.000  0.000   0.000  0.000   0.000  0.000   0.000  0.000   0.329
charge 0.769          0.000          0.000          0.000          0.000          0.000          0.000          0.000          0.000          0.769
            "#;

        let pdos = Procar::parse_projections_aux(txt, 4, false)?;
        assert_eq!(pdos.shape(), &[4, 9]);
        assert_eq!(pdos.column(0), arr1(&[0.060, 0.059, 0.339, 0.339]));
        assert_eq!(pdos.column(8), arr1(&[0.000, 0.000, 0.000, 0.000]));

        let txt: &str = r#"
band     1 # energy  -16.82688427 # occ.  1.00000000

ion      s     py     pz     px    dxy    dyz    dz2    dxz  x2-y2    tot
    1  0.306  0.000  0.003  0.000  0.000  0.000  0.000  0.000  0.000  0.309
    2  0.306  0.000  0.003  0.000  0.000  0.000  0.000  0.000  0.000  0.310
    3  0.112  0.000  0.000  0.000  0.000  0.000  0.000  0.000  0.000  0.112
tot    0.724  0.000  0.007  0.000  0.000  0.000  0.000  0.000  0.000  0.731
    1  0.013  0.000  0.000 -0.000 -0.000 -0.000  0.000  0.000  0.000  0.014
    2  0.013  0.000  0.000 -0.000 -0.000 -0.000  0.000  0.000  0.000  0.014
    3  0.005  0.000 -0.000  0.000  0.000  0.000  0.000 -0.000 -0.000  0.005
tot    0.032  0.000  0.000 -0.000  0.000  0.000  0.000 -0.000 -0.000  0.032
    1 -0.002  0.000 -0.000 -0.000  0.000 -0.000 -0.000  0.000  0.000 -0.002
    2 -0.002  0.000 -0.000 -0.000  0.000 -0.000 -0.000  0.000  0.000 -0.002
    3 -0.001 -0.000 -0.000  0.000 -0.000 -0.000 -0.000 -0.000 -0.000 -0.001
tot   -0.006  0.000 -0.000 -0.000 -0.000 -0.000 -0.000 -0.000 -0.000 -0.006
    1 -0.306  0.000 -0.003  0.000  0.000  0.000 -0.000  0.000  0.000 -0.309
    2 -0.306  0.000 -0.003  0.000  0.000  0.000 -0.000  0.000  0.000 -0.309
    3 -0.112  0.000 -0.000  0.000  0.000  0.000 -0.000  0.000  0.000 -0.112
tot   -0.723  0.000 -0.007  0.000  0.000  0.000 -0.000  0.000  0.000 -0.730
            "#;

        let pdos = Procar::parse_projections_aux(txt, 3, true)?;
        assert_eq!(pdos.shape(), &[12, 9]);
        assert_eq!(pdos.column(0), arr1(&[ 0.306,  0.306,  0.112,
                                           0.013,  0.013,  0.005,
                                          -0.002, -0.002, -0.001,
                                          -0.306, -0.306, -0.112]));
        assert_eq!(pdos.column(8), arr1(&[0.000, 0.000, 0.000,
                                          0.000, 0.000, 0.000,
                                          0.000, 0.000, 0.000,
                                          0.000, 0.000, 0.000]));
        Ok(())
    }

    #[test]
    fn test_parse_projected() -> Result<()> {
        let pdos = Procar::parse_projections(TXT_SMP, 2, 10, 10, 3, 3, false)?;
        assert_eq!(pdos.shape(), &[2, 10, 10, 3, 3]);
        assert_eq!(pdos.slice(s![0, 0, 0, .., ..]), arr2(&[[0.025, 0.002, 0.000],
                                                           [0.025, 0.002, 0.000],
                                                           [0.741, 0.000, 0.000]]));
        assert_eq!(pdos.slice(s![-1, -1, -1, .., ..]), arr2(&[[0.120, 0.153, 0.000],
                                                              [0.120, 0.153, 0.000],
                                                              [0.009, 0.052, 0.000]]));

        let pdos = Procar::parse_projections(TXT_NCL, 1, 81, 28, 3, 9, true)?;
        assert_eq!(pdos.shape(), &[4, 81, 28, 3, 9]);
        assert_eq!(pdos.slice(s![..,  0,  0,  0, 0]), arr1(&[0.306, 0.013, -0.002, -0.306]));
        assert_eq!(pdos.slice(s![.., -1, -1, -1, 0]), arr1(&[0.012, 0.000, -0.000, -0.012]));

        let pdos = Procar::parse_projections(TXT_STD, 2, 18, 49, 8, 9, false)?;
        assert_eq!(pdos.shape(), &[2, 18, 49, 8, 9]);
        assert_eq!(pdos.slice(s![ 0,  0,  0, .., 1]), arr1(&[0.498, 0.000, 0.498, 0.000, 0.000, 0.000, 0.000, 0.000]));
        assert_eq!(pdos.slice(s![-1, -1, -1, .., 0]), arr1(&[0.003, 0.000, 0.003, 0.000, 0.033, 0.000, 0.030, 0.000]));

        let pdos = Procar::parse_projections(TXT_PHS, 2, 60, 12, 3, 9, false)?;
        assert_eq!(pdos.shape(), &[2, 60, 12, 3, 9]);
        assert_eq!(pdos.slice(s![ 0,  0,  0, .., 0]), arr1(&[0.496, 0.496, 0.001]));
        assert_eq!(pdos.slice(s![-1, -1, -1, .., 3]), arr1(&[0.051, 0.051, 0.023]));

        let pdos = Procar::parse_projections(TXT_544, 2, 13, 20, 4, 9, false)?;
        assert_eq!(pdos.shape(), &[2, 13, 20, 4, 9]);
        assert_eq!(pdos.slice(s![ 0,  0,  0, .., 0]), arr1(&[0.060, 0.059, 0.339, 0.339]));
        assert_eq!(pdos.slice(s![-1, -1, -1, .., 0]), arr1(&[0.218, 0.000, 0.079, 0.079]));
        Ok(())
    }

    #[test]
    fn test_procar_from_file() -> Result<()> {
        Procar::from_file("tests/PROCAR")?;
        Procar::from_file("tests/PROCAR.simple")?;
        Procar::from_file("tests/PROCAR.phase")?;
        Procar::from_file("tests/PROCAR.ncl")?;
        Procar::from_file("tests/PROCAR.new_format_5.4.4")?;

        Ok(())
    }

}
