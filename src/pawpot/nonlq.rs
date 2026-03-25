// Reciprocal-space nonlocal projector (nonlq).

use std::f64::consts::TAU;
use num::complex::Complex;

use super::sph_harm::sph_r;
use super::gvec::{gvectors, gvectors_gamma};
use super::pawpotcar::PawPotcar;
use super::poscar::Poscar;

type C128 = Complex<f64>;

/// Reciprocal-space nonlocal projector.
pub struct Nonlq {
    /// Total number of plane waves.
    pub nplw: usize,
    /// G-vectors in fractional coordinates, shape (nplw,).
    pub gvec: Vec<[i32; 3]>,
    /// (G+k) in Cartesian (Å^{-1}), shape (nplw,).
    pub gk: Vec<[f64; 3]>,
    /// |G+k| in Å^{-1}, shape (nplw,).
    pub glen: Vec<f64>,
    /// Phase factor exp(i * 2π * G · R_atom) for each G and each atom.
    crexp: Vec<Vec<C128>>,
    /// Reciprocal projector * i^l for each element type.
    qproj: Vec<Vec<Vec<C128>>>,
    /// Element type index per atom.
    element_idx: Vec<usize>,
    /// Total lmmax per element type.
    lmmax_per_type: Vec<usize>,
    /// Total lmmax (sum over all atoms).
    pub total_lmmax: usize,
    /// Whether this is a gamma-only calculation.
    lgam: bool,
}

const TPI: f64 = TAU;

impl Nonlq {
    /// Build the nonlq projector.
    pub fn new(
        poscar: &Poscar,
        pawpot: &PawPotcar,
        encut: f64,
        kvec: &[f64; 3],
        bcell: &[[f64; 3]; 3],
        lgam: bool,
    ) -> Self {
        let ngrid = compute_ngrid(&poscar.cell, encut);
        let (gvec, gk, glen) = if lgam {
            gvectors_gamma(bcell, encut, kvec, &ngrid)
        } else {
            gvectors(bcell, encut, kvec, &ngrid)
        };
        let nplw = gvec.len();

        let natoms = poscar.natoms;
        let element_idx = poscar.element_idx();

        let scaled_positions = &poscar.positions;
        let mut crexp = vec![vec![C128::new(0.0, 0.0); natoms]; nplw];
        for (ig, g) in gvec.iter().enumerate() {
            for (iatom, pos) in scaled_positions.iter().enumerate() {
                let dot = g[0] as f64 * pos[0] + g[1] as f64 * pos[1] + g[2] as f64 * pos[2];
                crexp[ig][iatom] = Complex::new(0.0, TPI * dot).exp();
            }
        }

        let volume = poscar.volume().abs();

        let ntypes = pawpot.len();
        let mut qproj: Vec<Vec<Vec<C128>>> = Vec::with_capacity(ntypes);
        let mut lmmax_per_type = Vec::with_capacity(ntypes);

        for itype in 0..ntypes {
            let pp = pawpot.get(itype).unwrap();
            let lmmax = pp.lmmax();
            lmmax_per_type.push(lmmax);
            let gmax = pp.gmax();

            let lmax_pp = pp.ls().iter().copied().max().unwrap_or(0);
            let ylm_per_l: Vec<Vec<Vec<f64>>> = (0..=lmax_pp)
                .map(|l| sph_r(&gk, l))
                .collect();

            let mut tmp: Vec<Vec<C128>> = vec![vec![C128::new(0.0, 0.0); nplw]; lmmax];

            let spl_q = pp.spl_qproj();
            let ls = pp.ls();

            let mut ilm = 0usize;
            for (&l, spl) in ls.iter().zip(spl_q.iter()) {
                let tlp1 = 2 * l + 1;
                let il_factor = C128::new(0.0, 1.0).powi(l as i32);

                for (ig, &glen_ig) in glen.iter().enumerate().take(nplw) {
                    if glen_ig > gmax {
                        continue;
                    }
                    let q_val = spl.eval(glen_ig).unwrap_or(0.0);
                    let ylm_vals = &ylm_per_l[l][ig];

                    for (im, &y) in ylm_vals.iter().enumerate() {
                        let c = C128::new(q_val * y, 0.0) * il_factor;
                        tmp[ilm + im][ig] = c;
                    }
                }
                ilm += tlp1;
            }

            let inv_sqrt_vol = 1.0 / volume.sqrt();
            for row in tmp.iter_mut() {
                for c in row.iter_mut() {
                    *c *= inv_sqrt_vol;
                }
            }

            qproj.push(tmp);
        }

        let total_lmmax: usize = element_idx.iter().map(|&t| lmmax_per_type[t]).sum();

        if lgam {
            let sqrt2 = std::f64::consts::SQRT_2;
            for qp in qproj.iter_mut() {
                for row in qp.iter_mut() {
                    for c in row[1..].iter_mut() {
                        *c *= sqrt2;
                    }
                }
            }
        }

        Nonlq {
            nplw,
            gvec,
            gk,
            glen,
            crexp,
            qproj,
            element_idx,
            lmmax_per_type,
            total_lmmax,
            lgam,
        }
    }

    /// Project plane-wave coefficients onto all (lm, atom) channels.
    pub fn proj(&self, cptwf: &[C128]) -> Vec<C128> {
        assert_eq!(cptwf.len(), self.nplw, "nplw mismatch");
        let natoms = self.element_idx.len();
        let mut beta = Vec::with_capacity(self.total_lmmax);

        for iatom in 0..natoms {
            let ntype = self.element_idx[iatom];
            let lmmax = self.lmmax_per_type[ntype];
            let qp = &self.qproj[ntype];

            for qp_row in qp.iter().take(lmmax) {
                let mut sum = C128::new(0.0, 0.0);
                for ig in 0..self.nplw {
                    sum += cptwf[ig] * self.crexp[ig][iatom] * qp_row[ig];
                }
                let b = if self.lgam { C128::new(sum.re, 0.0) } else { sum };
                beta.push(b);
            }
        }
        beta
    }
}

fn compute_ngrid(cell: &[[f64; 3]; 3], encut: f64) -> [usize; 3] {
    const AUTOA: f64 = 0.529177249;
    const RYTOEV: f64 = 13.605826;
    let anorm = [
        vec_norm(&cell[0]),
        vec_norm(&cell[1]),
        vec_norm(&cell[2]),
    ];
    let cutof: [f64; 3] = std::array::from_fn(|i| {
        ((encut / RYTOEV).sqrt() / (TPI / (anorm[i] / AUTOA))).ceil()
    });
    [
        (2.0 * cutof[0] + 1.0) as usize,
        (2.0 * cutof[1] + 1.0) as usize,
        (2.0 * cutof[2] + 1.0) as usize,
    ]
}

fn vec_norm(v: &[f64; 3]) -> f64 {
    (v[0] * v[0] + v[1] * v[1] + v[2] * v[2]).sqrt()
}
