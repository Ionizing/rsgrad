use std::fs;
use std::path::Path;

use anyhow::{
    Context,
    Result,
    bail,
};

use super::cspline::{CubicSpline, BoundaryCondition};


#[derive(Clone)]
pub struct PawPotcarSingle {
    symbol:         String,
    zval:           f64,
    rmax:           f64,
    gmax:           f64,
    radngrid:       usize,
    ndata:          usize,
    ls:             Vec<usize>,
    radgrid:        Vec<f64>,
    proj_rgrid:     Vec<f64>,
    proj_qgrid:     Vec<f64>,
    aepotential:    Vec<f64>,
    aecorecharge:   Vec<f64>,
    pswaves:        Vec<Vec<f64>>,
    aewaves:        Vec<Vec<f64>>,
    reciprojs:      Vec<Vec<f64>>,
    realprojs:      Vec<Vec<f64>>,
    simp_weights:   Vec<f64>,
    spl_qproj:      Vec<CubicSpline>,
    spl_rproj:      Vec<CubicSpline>,
}


impl std::str::FromStr for PawPotcarSingle {
    type Err = anyhow::Error;
    fn from_str(s: &str) -> Result<Self> {
        Self::parse(s)
    }
}

impl PawPotcarSingle {
    pub fn from_file<T>(fname: &T) -> Result<Self>
    where T: AsRef<Path> + ?Sized {
        let s = fs::read_to_string(fname)?;
        s.parse::<Self>()
    }

    fn parse(s: &str) -> Result<Self> {
        if s.contains("AE potential") {
            bail!("This POTCAR is not supported: AE potential.");
        }

        let symbol       = Self::parse_symbol(s)?;
        let zval         = Self::parse_zval(s)?;
        let (ls, rmax)   = Self::parse_ls_rmax(s)?;
        let gmax         = Self::parse_gmax(s)?;
        let radngrid     = Self::parse_radial_ngrid(s)?;
        let ndata        = Self::parse_ndata(s)?;
        let radgrid      = Self::parse_radgrid(s, radngrid)?;
        let proj_rgrid   = Self::gen_projgrid(ndata, rmax);
        let proj_qgrid   = Self::gen_projgrid(ndata, gmax);
        let aepotential  = Self::parse_aepotential(s, radngrid)?;
        let aecorecharge = Self::parse_aecorecharge(s, radngrid)?;
        let pswaves      = Self::parse_pswaves(s, ls.len(), radngrid)?;
        let aewaves      = Self::parse_aewaves(s, ls.len(), radngrid)?;
        let realprojs    = Self::parse_realprojs(s, ls.len(), ndata)?;
        let reciprojs    = Self::parse_reciprojs(s, ls.len(), ndata)?;
        let simp_weights = Self::gen_simp_weights(&radgrid);
        let spl_qproj    = Self::gen_qsplines(&proj_qgrid, &reciprojs);
        let spl_rproj    = Self::gen_rsplines(&ls, &proj_rgrid, &realprojs, rmax, ndata);

        Ok(Self {
            symbol,
            zval,
            rmax,
            gmax,
            radngrid,
            ndata,
            ls,
            radgrid,
            proj_rgrid,
            proj_qgrid,
            aepotential,
            aecorecharge,
            pswaves,
            aewaves,
            reciprojs,
            realprojs,
            simp_weights,
            spl_qproj,
            spl_rproj,
        })
    }

    // -------------------------------------------------------------------------
    // Public accessors
    // -------------------------------------------------------------------------

    pub fn symbol(&self) -> &str { &self.symbol }
    pub fn zval(&self) -> f64 { self.zval }
    pub fn rmax(&self) -> f64 { self.rmax }
    pub fn gmax(&self) -> f64 { self.gmax }
    pub fn radngrid(&self) -> usize { self.radngrid }
    pub fn ndata(&self) -> usize { self.ndata }
    pub fn ls(&self) -> &[usize] { &self.ls }
    pub fn radgrid(&self) -> &[f64] { &self.radgrid }
    pub fn proj_rgrid(&self) -> &[f64] { &self.proj_rgrid }
    pub fn proj_qgrid(&self) -> &[f64] { &self.proj_qgrid }
    pub fn aepotential(&self) -> &[f64] { &self.aepotential }
    pub fn aecorecharge(&self) -> &[f64] { &self.aecorecharge }
    pub fn pswaves(&self) -> &[Vec<f64>] { &self.pswaves }
    pub fn aewaves(&self) -> &[Vec<f64>] { &self.aewaves }
    pub fn realprojs(&self) -> &[Vec<f64>] { &self.realprojs }
    pub fn reciprojs(&self) -> &[Vec<f64>] { &self.reciprojs }
    pub fn simp_weights(&self) -> &[f64] { &self.simp_weights }
    pub fn spl_qproj(&self) -> &[CubicSpline] { &self.spl_qproj }
    pub fn spl_rproj(&self) -> &[CubicSpline] { &self.spl_rproj }

    // -------------------------------------------------------------------------
    // Computed properties
    // -------------------------------------------------------------------------

    pub fn lmax(&self) -> usize { self.ls.len() }

    pub fn lmmax(&self) -> usize {
        self.ls.iter().map(|&l| 2 * l + 1).sum()
    }

    pub fn ilm(&self) -> Vec<(usize, usize, i64)> {
        self.ls.iter().enumerate()
            .flat_map(|(i, &l)| {
                let li = l as i64;
                (-li ..= li).map(move |m| (i, l, m))
            })
            .collect()
    }

    // -------------------------------------------------------------------------
    // Integration
    // -------------------------------------------------------------------------

    pub fn radial_simp_integrate(&self, f: &[f64]) -> f64 {
        assert_eq!(self.radngrid, f.len(),
            "radial_simp_integrate: f length {} != radngrid {}", f.len(), self.radngrid);
        self.simp_weights.iter()
            .zip(f.iter())
            .map(|(w, y)| w * y)
            .sum()
    }

    // -------------------------------------------------------------------------
    // PAW correction matrices
    // -------------------------------------------------------------------------

    pub fn get_qij(&self) -> Vec<Vec<f64>> {
        let lmmax = self.lmmax();
        let ilm   = self.ilm();
        let mut qij = vec![vec![0.0f64; lmmax]; lmmax];

        for ii in 0 .. lmmax {
            for jj in 0 ..= ii {
                let (n1, l1, m1) = ilm[ii];
                let (n2, l2, m2) = ilm[jj];
                if l1 == l2 && m1 == m2 {
                    let integrand: Vec<f64> = self.aewaves[n1].iter()
                        .zip(self.aewaves[n2].iter())
                        .zip(self.pswaves[n1].iter())
                        .zip(self.pswaves[n2].iter())
                        .map(|(((ae1, ae2), ps1), ps2)| ae1 * ae2 - ps1 * ps2)
                        .collect();
                    let val = self.radial_simp_integrate(&integrand);
                    qij[ii][jj] = val;
                    qij[jj][ii] = val;
                }
            }
        }

        qij
    }

    pub fn get_nablaij(&self) -> Vec<Vec<Vec<f64>>> {
        let lmmax = self.lmmax();
        let ilm   = self.ilm();
        let rr    = &self.radgrid;

        let grad_ae: Vec<Vec<f64>> = self.aewaves.iter().map(|w| radial_grad(rr, w)).collect();
        let grad_ps: Vec<Vec<f64>> = self.pswaves.iter().map(|w| radial_grad(rr, w)).collect();

        let mut nablaij = vec![vec![vec![0.0f64; lmmax]; lmmax]; 3];

        for ii in 0..lmmax {
            for jj in 0..lmmax {
                let (n1, l1, m1) = ilm[ii];
                let (n2, l2, m2) = ilm[jj];
                let l1i = l1 as i64;
                let m1i = m1;
                let l2i = l2 as i64;
                let m2i = m2;

                let a1 = gaunt_a1(l1i, m1i, l2i, m2i);
                let a2 = ylm_nabla_rlylm_tab(l1i, m1i, l2i, m2i);

                if a1.iter().chain(a2.iter()).all(|&v| v.abs() < 1e-14) {
                    continue;
                }

                let l2p1 = (l2 + 1) as f64;
                let r1_ae = self.radial_simp_integrate(
                    &rr.iter().enumerate()
                        .map(|(k, &r)| self.aewaves[n1][k] * (grad_ae[n2][k] - l2p1 * self.aewaves[n2][k] / r))
                        .collect::<Vec<_>>(),
                );
                let r1_ps = self.radial_simp_integrate(
                    &rr.iter().enumerate()
                        .map(|(k, &r)| self.pswaves[n1][k] * (grad_ps[n2][k] - l2p1 * self.pswaves[n2][k] / r))
                        .collect::<Vec<_>>(),
                );
                let r1 = r1_ae - r1_ps;

                let r2_ae = self.radial_simp_integrate(
                    &rr.iter().enumerate()
                        .map(|(k, &r)| self.aewaves[n1][k] * self.aewaves[n2][k] / r)
                        .collect::<Vec<_>>(),
                );
                let r2_ps = self.radial_simp_integrate(
                    &rr.iter().enumerate()
                        .map(|(k, &r)| self.pswaves[n1][k] * self.pswaves[n2][k] / r)
                        .collect::<Vec<_>>(),
                );
                let r2 = r2_ae - r2_ps;

                for alpha in 0..3 {
                    nablaij[alpha][ii][jj] = r1 * a1[alpha] + r2 * a2[alpha];
                }
            }
        }

        nablaij
    }

    // -------------------------------------------------------------------------
    // Private helpers: spline construction
    // -------------------------------------------------------------------------

    fn gen_qsplines(proj_qgrid: &[f64], reciprojs: &[Vec<f64>]) -> Vec<CubicSpline> {
        reciprojs.iter()
            .map(|qproj| CubicSpline::from_xy_natural(proj_qgrid, qproj))
            .collect()
    }

    fn gen_rsplines(
        ls: &[usize],
        proj_rgrid: &[f64],
        realprojs: &[Vec<f64>],
        rmax: f64,
        ndata: usize,
    ) -> Vec<CubicSpline> {
        let h = rmax / ndata as f64;
        ls.iter().zip(realprojs.iter())
            .map(|(&l, rproj)| {
                let y1p = if l == 1 {
                    (rproj[1] - rproj[0]) / h
                } else {
                    0.0
                };
                let bc = [
                    BoundaryCondition::FirstDerivative(y1p),
                    BoundaryCondition::SecondDerivative(0.0),
                ];
                CubicSpline::from_xy_with_bc(proj_rgrid, rproj, &bc)
            })
            .collect()
    }

    // -------------------------------------------------------------------------
    // Private helpers: zval
    // -------------------------------------------------------------------------

    fn parse_zval(s: &str) -> Result<f64> {
        let zval = s.lines()
            .filter(|x| !x.trim().is_empty())
            .nth(1).with_context(|| "Cannot find zval in POTCAR.")?
            .trim()
            .parse::<f64>()?;
        Ok(zval)
    }

    fn parse_radgrid(s: &str, radngrid: usize) -> Result<Vec<f64>> {
        let radgrid = s.split("\n grid")
            .nth(1).with_context(|| "Cannot find grid data from POTCAR.")?
            .split_ascii_whitespace()
            .take(radngrid)
            .map(|x| x.parse::<f64>().unwrap())
            .collect::<Vec<_>>();
        Ok(radgrid)
    }

    fn parse_radial_ngrid(s: &str) -> Result<usize> {
        let ngrid = s.split("PAW radial sets")
            .nth(1).with_context(|| "Cannot find radial grid number.")?
            .split_ascii_whitespace()
            .nth(0).unwrap()
            .parse::<usize>()?;

        Ok(ngrid)
    }

    fn parse_symbol(s: &str) -> Result<String> {
        let symbol = s.lines().find(|x| !x.trim().is_empty()).unwrap()
            .split_ascii_whitespace()
            .nth(1).unwrap()
            .split('_')
            .next().unwrap()
            .to_owned();

        Ok(symbol)
    }

    fn parse_gmax(s: &str) -> Result<f64> {
        let gmax = s.split(" Non local Part")
            .next().unwrap()
            .lines()
            .last().unwrap()
            .split_ascii_whitespace()
            .next().unwrap()
            .parse::<f64>().unwrap();
        Ok(gmax)
    }

    fn parse_ls_rmax(s: &str) -> Result<(Vec<usize>, f64)> {
        let mut lss = Vec::<usize>::with_capacity(4);
        let mut rmax = 0.0f64;

        s.split(" Non local Part")
            .skip(1)
            .for_each(|x| {
                let mut it = x.split_ascii_whitespace();
                let l      = it.next().unwrap().parse::<usize>().unwrap();
                let n      = it.next().unwrap().parse::<usize>().unwrap();
                let _rmax  = it.next().unwrap().parse::<f64>().unwrap();
                lss.extend(vec![l; n]);
                rmax = _rmax;
            });

        Ok((lss, rmax))
    }

    fn parse_ndata(s: &str) -> Result<usize> {
        let ndata = s.split("   NDATA  =")
            .nth(1).with_context(|| "Cannot find NDATA info.")?
            .split_ascii_whitespace()
            .next().unwrap()
            .parse::<usize>()?;

        Ok(ndata)
    }

    fn gen_projgrid(ndata: usize, rmax: f64) -> Vec<f64> {
        assert!(rmax > 0.0, "PawPot.gen_projgrid: rmax cannot be less than 0.0.");
        (0 .. ndata)
            .map(|x| x as f64 * rmax / ndata as f64)
            .collect()
    }

    fn parse_aepotential(s: &str, radngrid: usize) -> Result<Vec<f64>> {
        let aepotential = s.split("\n aepotential")
            .nth(1).with_context(|| "Cannot find aepotential data.")?
            .split_ascii_whitespace()
            .take(radngrid)
            .map(|x| x.parse::<f64>().unwrap())
            .collect::<Vec<_>>();

        Ok(aepotential)
    }

    fn parse_aecorecharge(s: &str, radngrid: usize) -> Result<Vec<f64>> {
        let aecorecharge = s.split("\n core charge-density\n")
            .nth(1).with_context(|| "Cannot find core charge density data.")?
            .split_ascii_whitespace()
            .take(radngrid)
            .map(|x| x.parse::<f64>().unwrap())
            .collect::<Vec<_>>();

        Ok(aecorecharge)
    }

    fn parse_pswaves(s: &str, nls: usize, radngrid: usize) -> Result<Vec<Vec<f64>>> {
        let pswaves = s.split("pseudo wavefunction")
            .skip(1)
            .map(|xs| {
                xs.split_ascii_whitespace()
                    .take(radngrid)
                    .map(|x| x.parse::<f64>().unwrap())
                    .collect::<Vec<_>>()
            })
            .collect::<Vec<_>>();

        assert_eq!(pswaves.len(), nls, "Pseudo wavefunction sets not complete.");
        Ok(pswaves)
    }

    fn parse_aewaves(s: &str, nls: usize, radngrid: usize) -> Result<Vec<Vec<f64>>> {
        let aewaves = s.split("ae wavefunction")
            .skip(1)
            .map(|xs| {
                xs.split_ascii_whitespace()
                    .take(radngrid)
                    .map(|x| x.parse::<f64>().unwrap())
                    .collect::<Vec<_>>()
            })
            .collect::<Vec<_>>();

        assert_eq!(aewaves.len(), nls, "Ae wavefunction sets not complete.");
        Ok(aewaves)
    }

    fn parse_realprojs(s: &str, nls: usize, ndata: usize) -> Result<Vec<Vec<f64>>> {
        let realprojs = s.split(" Real Space Part")
            .skip(1)
            .map(|xs| {
                xs.split_ascii_whitespace()
                    .take(ndata)
                    .map(|x| x.parse::<f64>().unwrap())
                    .collect::<Vec<_>>()
            })
            .collect::<Vec<_>>();

        assert_eq!(realprojs.len(), nls, "Realprojs sets not complete.");
        Ok(realprojs)
    }

    fn parse_reciprojs(s: &str, nls: usize, ndata: usize) -> Result<Vec<Vec<f64>>> {
        let reciprojs = s.split(" Reciprocal Space Part")
            .skip(1)
            .map(|xs| {
                xs.split_ascii_whitespace()
                    .take(ndata)
                    .map(|x| x.parse::<f64>().unwrap())
                    .collect::<Vec<_>>()
            })
            .collect::<Vec<_>>();

        assert_eq!(reciprojs.len(), nls, "Reciprojs sets not complete.");
        Ok(reciprojs)
    }

    fn gen_simp_weights(radgrid: &[f64]) -> Vec<f64> {
        let     ngrid      = radgrid.len();
        let mut rad_simp_w = vec![0.0f64; ngrid];
        let     h          = f64::ln(
            (radgrid[ngrid-1] / radgrid[0]).powf(1.0 / (ngrid - 1) as f64)
            );

        for i in (3 - (ngrid % 2) ..= ngrid-1).rev().step_by(2) {
            rad_simp_w[i  ] += h * radgrid[i  ] / 3.0;
            rad_simp_w[i-1] = h * radgrid[i-1] * 4.0 / 3.0;
            rad_simp_w[i-2] = h * radgrid[i-2] / 3.0;
        }

        rad_simp_w
    }
}



pub struct PawPotcar {
    pots: Vec<PawPotcarSingle>,
}


impl std::str::FromStr for PawPotcar {
    type Err = anyhow::Error;
    fn from_str(s: &str) -> Result<Self> {
        let pots = s.split(" End of Dataset")
            .filter(|x| !x.trim().is_empty())
            .map(|x: &str| x.parse::<PawPotcarSingle>())
            .collect::<Result<Vec<PawPotcarSingle>>>()?;

        Ok(Self { pots })
    }
}

impl PawPotcar {
    pub fn from_file<T>(fname: &T) -> Result<Self>
    where T: AsRef<Path> + ?Sized {
        let s = fs::read_to_string(fname)?;
        s.parse::<Self>()
    }

    pub fn len(&self) -> usize { self.pots.len() }

    pub fn is_empty(&self) -> bool { self.pots.is_empty() }

    pub fn get(&self, i: usize) -> Option<&PawPotcarSingle> { self.pots.get(i) }

    pub fn get_by_symbol(&self, sym: &str) -> Option<&PawPotcarSingle> {
        self.pots.iter().find(|p| p.symbol() == sym)
    }

    pub fn iter(&self) -> std::slice::Iter<'_, PawPotcarSingle> { self.pots.iter() }
}

impl std::ops::Index<usize> for PawPotcar {
    type Output = PawPotcarSingle;
    fn index(&self, i: usize) -> &Self::Output { &self.pots[i] }
}

// ---------------------------------------------------------------------------
// Nabla matrix element helpers
// ---------------------------------------------------------------------------

fn radial_grad(rr: &[f64], ff: &[f64]) -> Vec<f64> {
    let n = rr.len();
    assert!(n >= 2);
    let mut g = vec![0.0f64; n];
    g[0] = (ff[1] - ff[0]) / (rr[1] - rr[0]);
    for i in 1..n - 1 {
        g[i] = (ff[i + 1] - ff[i - 1]) / (rr[i + 1] - rr[i - 1]);
    }
    g[n - 1] = (ff[n - 1] - ff[n - 2]) / (rr[n - 1] - rr[n - 2]);
    g
}

#[allow(clippy::excessive_precision)]
fn gaunt_a1(l1: i64, m1: i64, l2: i64, m2: i64) -> [f64; 3] {
    match (l1, m1, l2, m2) {
        (0, 0, 1,-1) => [0.0, 0.577350269189625, 0.0],
        (0, 0, 1, 0) => [0.0, 0.0, 0.577350269189626],
        (0, 0, 1, 1) => [0.577350269189625, 0.0, 0.0],
        (1,-1, 0, 0) => [0.0, 0.577350269189625, 0.0],
        (1, 0, 0, 0) => [0.0, 0.0, 0.577350269189626],
        (1, 1, 0, 0) => [0.577350269189625, 0.0, 0.0],
        (1,-1, 2,-2) => [0.447213595499958, 0.0, 0.0],
        (1,-1, 2,-1) => [0.0, 0.0, 0.447213595499958],
        (1,-1, 2, 0) => [0.0,-0.258198889747161, 0.0],
        (1,-1, 2, 2) => [0.0,-0.447213595499958, 0.0],
        (1, 0, 2,-1) => [0.0, 0.447213595499958, 0.0],
        (1, 0, 2, 0) => [0.0, 0.0, 0.516397779494322],
        (1, 0, 2, 1) => [0.447213595499958, 0.0, 0.0],
        (1, 1, 2,-2) => [0.0, 0.447213595499958, 0.0],
        (1, 1, 2, 0) => [-0.258198889747161, 0.0, 0.0],
        (1, 1, 2, 1) => [0.0, 0.0, 0.447213595499958],
        (1, 1, 2, 2) => [0.447213595499958, 0.0, 0.0],
        (2,-2, 1,-1) => [0.447213595499958, 0.0, 0.0],
        (2,-2, 1, 1) => [0.0, 0.447213595499958, 0.0],
        (2,-1, 1,-1) => [0.0, 0.0, 0.447213595499958],
        (2,-1, 1, 0) => [0.0, 0.447213595499958, 0.0],
        (2, 0, 1,-1) => [0.0,-0.258198889747161, 0.0],
        (2, 0, 1, 0) => [0.0, 0.0, 0.516397779494322],
        (2, 0, 1, 1) => [-0.258198889747161, 0.0, 0.0],
        (2, 1, 1, 0) => [0.447213595499958, 0.0, 0.0],
        (2, 1, 1, 1) => [0.0, 0.0, 0.447213595499958],
        (2, 2, 1,-1) => [0.0,-0.447213595499958, 0.0],
        (2, 2, 1, 1) => [0.447213595499958, 0.0, 0.0],
        (2,-2, 3,-3) => [0.462910049886276, 0.0, 0.0],
        (2,-2, 3,-2) => [0.0, 0.0, 0.377964473009227],
        (2,-2, 3,-1) => [-0.119522860933439, 0.0, 0.0],
        (2,-2, 3, 1) => [0.0,-0.119522860933439, 0.0],
        (2,-2, 3, 3) => [0.0,-0.462910049886276, 0.0],
        (2,-1, 3,-2) => [0.377964473009227, 0.0, 0.0],
        (2,-1, 3,-1) => [0.0, 0.0, 0.478091443733757],
        (2,-1, 3, 0) => [0.0,-0.292770021884560, 0.0],
        (2,-1, 3, 2) => [0.0,-0.377964473009227, 0.0],
        (2, 0, 3,-1) => [0.0, 0.414039335605412, 0.0],
        (2, 0, 3, 0) => [0.0, 0.0, 0.507092552837110],
        (2, 0, 3, 1) => [0.414039335605412, 0.0, 0.0],
        (2, 1, 3,-2) => [0.0, 0.377964473009227, 0.0],
        (2, 1, 3, 0) => [-0.292770021884560, 0.0, 0.0],
        (2, 1, 3, 1) => [0.0, 0.0, 0.478091443733757],
        (2, 1, 3, 2) => [0.377964473009227, 0.0, 0.0],
        (2, 2, 3,-3) => [0.0, 0.462910049886276, 0.0],
        (2, 2, 3,-1) => [0.0, 0.119522860933439, 0.0],
        (2, 2, 3, 1) => [-0.119522860933439, 0.0, 0.0],
        (2, 2, 3, 2) => [0.0, 0.0, 0.377964473009227],
        (2, 2, 3, 3) => [0.462910049886276, 0.0, 0.0],
        (3,-3, 2,-2) => [0.462910049886276, 0.0, 0.0],
        (3,-3, 2, 2) => [0.0, 0.462910049886276, 0.0],
        (3,-2, 2,-2) => [0.0, 0.0, 0.377964473009227],
        (3,-2, 2,-1) => [0.377964473009227, 0.0, 0.0],
        (3,-2, 2, 1) => [0.0, 0.377964473009227, 0.0],
        (3,-1, 2,-2) => [-0.119522860933439, 0.0, 0.0],
        (3,-1, 2,-1) => [0.0, 0.0, 0.478091443733757],
        (3,-1, 2, 0) => [0.0, 0.414039335605412, 0.0],
        (3,-1, 2, 2) => [0.0, 0.119522860933439, 0.0],
        (3, 0, 2,-1) => [0.0,-0.292770021884560, 0.0],
        (3, 0, 2, 0) => [0.0, 0.0, 0.507092552837110],
        (3, 0, 2, 1) => [-0.292770021884560, 0.0, 0.0],
        (3, 1, 2,-2) => [0.0,-0.119522860933439, 0.0],
        (3, 1, 2, 0) => [0.414039335605412, 0.0, 0.0],
        (3, 1, 2, 1) => [0.0, 0.0, 0.478091443733757],
        (3, 1, 2, 2) => [-0.119522860933439, 0.0, 0.0],
        (3, 2, 2,-1) => [0.0,-0.377964473009227, 0.0],
        (3, 2, 2, 1) => [0.377964473009227, 0.0, 0.0],
        (3, 2, 2, 2) => [0.0, 0.0, 0.377964473009227],
        (3, 3, 2,-2) => [0.0,-0.462910049886276, 0.0],
        (3, 3, 2, 2) => [0.462910049886276, 0.0, 0.0],
        _ => [0.0, 0.0, 0.0],
    }
}

#[allow(clippy::excessive_precision)]
fn ylm_nabla_rlylm_tab(l1: i64, m1: i64, l2: i64, m2: i64) -> [f64; 3] {
    match (l1, m1, l2, m2) {
        (0, 0, 1,-1) => [0.0, 1.732050807568877, 0.0],
        (0, 0, 1, 0) => [0.0, 0.0, 1.732050807568877],
        (0, 0, 1, 1) => [1.732050807568877, 0.0, 0.0],
        (1,-1, 2,-2) => [2.236_067_977_499_79, 0.0, 0.0],
        (1,-1, 2,-1) => [0.0, 0.0, 2.236_067_977_499_79],
        (1,-1, 2, 0) => [0.0,-1.290994448735806, 0.0],
        (1,-1, 2, 2) => [0.0,-2.236_067_977_499_79, 0.0],
        (1, 0, 2,-1) => [0.0, 2.236_067_977_499_79, 0.0],
        (1, 0, 2, 0) => [0.0, 0.0, 2.581988897471611],
        (1, 0, 2, 1) => [2.236_067_977_499_79, 0.0, 0.0],
        (1, 1, 2,-2) => [0.0, 2.236_067_977_499_79, 0.0],
        (1, 1, 2, 0) => [-1.290994448735806, 0.0, 0.0],
        (1, 1, 2, 1) => [0.0, 0.0, 2.236_067_977_499_79],
        (1, 1, 2, 2) => [2.236_067_977_499_79, 0.0, 0.0],
        (2,-2, 3,-3) => [3.240_370_349_203_93, 0.0, 0.0],
        (2,-2, 3,-2) => [0.0, 0.0, 2.645_751_311_064_59],
        (2,-2, 3,-1) => [-0.836660026534076, 0.0, 0.0],
        (2,-2, 3, 1) => [0.0,-0.836660026534076, 0.0],
        (2,-2, 3, 3) => [0.0,-3.240_370_349_203_93, 0.0],
        (2,-1, 3,-2) => [2.645_751_311_064_59, 0.0, 0.0],
        (2,-1, 3,-1) => [0.0, 0.0, 3.346640106136303],
        (2,-1, 3, 0) => [0.0,-2.049_390_153_191_92, 0.0],
        (2,-1, 3, 2) => [0.0,-2.645_751_311_064_59, 0.0],
        (2, 0, 3,-1) => [0.0, 2.898275349237888, 0.0],
        (2, 0, 3, 0) => [0.0, 0.0, 3.549_647_869_859_77],
        (2, 0, 3, 1) => [2.898275349237888, 0.0, 0.0],
        (2, 1, 3,-2) => [0.0, 2.645_751_311_064_59, 0.0],
        (2, 1, 3, 0) => [-2.049_390_153_191_92, 0.0, 0.0],
        (2, 1, 3, 1) => [0.0, 0.0, 3.346640106136303],
        (2, 1, 3, 2) => [2.645_751_311_064_59, 0.0, 0.0],
        (2, 2, 3,-3) => [0.0, 3.240_370_349_203_93, 0.0],
        (2, 2, 3,-1) => [0.0, 0.836660026534076, 0.0],
        (2, 2, 3, 1) => [-0.836660026534076, 0.0, 0.0],
        (2, 2, 3, 2) => [0.0, 0.0, 2.645751311064591],
        (2, 2, 3, 3) => [3.240_370_349_203_93, 0.0, 0.0],
        _ => [0.0, 0.0, 0.0],
    }
}
