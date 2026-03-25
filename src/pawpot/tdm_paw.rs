// Transition Dipole Moment (TDM) calculator using the p-r relation with PAW correction.

use std::path::Path;

use super::npy::NpzWriter;

use anyhow::Result;
use num::complex::Complex;

use super::gvec::{gvectors, gvectors_gamma};
use super::nonlq::Nonlq;
use super::pawpotcar::PawPotcar;
use super::poscar::Poscar;
use super::wavecar::Wavecar;

type C128 = Complex<f64>;

// Physical constants
const AUTOA: f64 = 0.529177249;
const RYTOEV: f64 = 13.605826;
const AUTDEBYE: f64 = 2.541746; // e·Bohr to Debye (same as VASP/VaspBandUnfolding)
const TPI: f64 = std::f64::consts::TAU;

/// One TDM entry: initial band i → final band j.
#[derive(Debug)]
pub struct TdmPawEntry {
    pub iband: usize,
    pub jband: usize,
    pub ei: f64,
    pub ej: f64,
    pub de: f64,
    /// Dipole matrix element |d| in Debye, (x, y, z).
    pub dipole: [C128; 3],
}

/// Compute the TDM between `ibands` and `jbands` at a given spin / k-point.
#[allow(clippy::too_many_arguments)]
pub fn compute_tdm_paw(
    poscar: &Poscar,
    pawpot: &PawPotcar,
    wavecar: &mut Wavecar,
    ispin: usize,
    ikpt: usize,
    ibands: &[usize],
    jbands: &[usize],
    lpseudo: bool,
) -> Result<Vec<TdmPawEntry>> {
    let encut = wavecar.header.encut;
    let bcell = wavecar.bcell();
    let kvec  = wavecar.kpoints[ikpt - 1].kvec;

    let ngrid = ngrid_from_cell(&poscar.cell, encut);
    let nplw_stored = wavecar.kpoints[ikpt - 1].nplw;
    let nplw_full   = gvectors(&bcell, encut, &kvec, &ngrid).0.len();
    let nplw_gamma  = gvectors_gamma(&bcell, encut, &kvec, &ngrid).0.len();
    let lgam = nplw_stored == nplw_gamma && nplw_stored != nplw_full;

    let nonlq = Nonlq::new(poscar, pawpot, encut, &kvec, &bcell, lgam);

    let (_, gk, _) = if lgam {
        gvectors_gamma(&bcell, encut, &kvec, &ngrid)
    } else {
        gvectors(&bcell, encut, &kvec, &ngrid)
    };
    let nplw = gk.len();

    let nablaij: Vec<Vec<Vec<Vec<f64>>>> = (0..pawpot.len())
        .map(|it| pawpot.get(it).unwrap().get_nablaij())
        .collect();

    let energies: Vec<f64> = wavecar.kpoints[ikpt - 1].eigenvalues.clone();

    let read_coeff = |wavecar: &mut Wavecar, iband: usize| -> Result<Vec<C128>> {
        wavecar.read_band_coeff(ispin, ikpt, iband)
    };

    let mut coeff_cache: std::collections::HashMap<usize, Vec<C128>> =
        std::collections::HashMap::new();
    for &b in ibands.iter().chain(jbands.iter()) {
        if let std::collections::hash_map::Entry::Vacant(e) = coeff_cache.entry(b) {
            let c = read_coeff(wavecar, b)?;
            e.insert(c);
        }
    }

    let mut entries = Vec::new();

    for &ib in ibands {
        let ei = energies[ib - 1];
        let ci = &coeff_cache[&ib];

        let beta_i = if !lpseudo {
            Some(nonlq.proj(ci))
        } else {
            None
        };

        for &jb in jbands {
            let ej = energies[jb - 1];
            let de = ej - ei;
            let cj = &coeff_cache[&jb];

            let beta_j = if !lpseudo {
                Some(nonlq.proj(cj))
            } else {
                None
            };

            // PS momentum matrix: dp[α] = Σ_G conj(C_j,G) * C_i,G * (G+k)_α
            let mut dp = [C128::new(0.0, 0.0); 3];
            if lgam {
                for ig in 0..nplw {
                    let cjc = cj[ig].conj();
                    let cci = ci[ig].conj();
                    for (alpha, dp_alpha) in dp.iter_mut().enumerate() {
                        let g = gk[ig][alpha];
                        *dp_alpha += (cjc * ci[ig] - cj[ig] * cci) * g / 2.0;
                    }
                }
            } else {
                for ig in 0..nplw {
                    let cjc = cj[ig].conj();
                    for (alpha, dp_alpha) in dp.iter_mut().enumerate() {
                        *dp_alpha += cjc * ci[ig] * gk[ig][alpha];
                    }
                }
            }

            // PAW one-centre correction
            if !lpseudo {
                let bi = beta_i.as_ref().unwrap();
                let bj = beta_j.as_ref().unwrap();
                let element_idx = poscar.element_idx();
                let lmmax_per_type: Vec<usize> = (0..pawpot.len())
                    .map(|it| pawpot.get(it).unwrap().lmmax())
                    .collect();

                let mut nproj = 0usize;
                for &ntype in &element_idx {
                    let lmmax = lmmax_per_type[ntype];
                    let nab = &nablaij[ntype];

                    for alpha in 0..3 {
                        let mut corr = C128::new(0.0, 0.0);
                        for m in 0..lmmax {
                            let mut tmp = C128::new(0.0, 0.0);
                            for n in 0..lmmax {
                                tmp += nab[alpha][m][n] * bi[nproj + n];
                            }
                            corr += bj[nproj + m].conj() * tmp;
                        }
                        // correction += -i * corr
                        dp[alpha] += C128::new(corr.im, -corr.re);
                    }
                    nproj += lmmax;
                }
            }

            // p-r relation: d[α] = -i / (dE/(2*RYTOEV)) * dp[α] * AUTOA * AUTDEBYE
            let dipole = if de.abs() < 1e-6 {
                [C128::new(0.0, 0.0); 3]
            } else {
                let fac = -1.0 / (de / (2.0 * RYTOEV)) * AUTOA * AUTDEBYE;
                std::array::from_fn(|alpha| C128::new(dp[alpha].im * fac, -dp[alpha].re * fac))
            };

            entries.push(TdmPawEntry { iband: ib, jband: jb, ei, ej, de, dipole });
        }
    }

    Ok(entries)
}

/// Write TDM results to a `.npz` file.
pub fn write_tdm_paw_npz(path: impl AsRef<Path>, entries: &[TdmPawEntry]) -> Result<()> {
    let n = entries.len();
    let iband: Vec<i64>  = entries.iter().map(|e| e.iband as i64).collect();
    let jband: Vec<i64>  = entries.iter().map(|e| e.jband as i64).collect();
    let ei:    Vec<f64>  = entries.iter().map(|e| e.ei).collect();
    let ej:    Vec<f64>  = entries.iter().map(|e| e.ej).collect();
    let de:    Vec<f64>  = entries.iter().map(|e| e.de).collect();
    let dipole: Vec<C128> = entries.iter()
        .flat_map(|e| e.dipole.iter().copied())
        .collect();

    let mut npz = NpzWriter::create(path.as_ref())?;
    npz.write_i64("iband",  &iband,  &[n as u64])?;
    npz.write_i64("jband",  &jband,  &[n as u64])?;
    npz.write_f64("E_i",    &ei,     &[n as u64])?;
    npz.write_f64("E_j",    &ej,     &[n as u64])?;
    npz.write_f64("dE",     &de,     &[n as u64])?;
    npz.write_c128("dipole", &dipole, &[n as u64, 3])?;
    npz.finish()?;
    Ok(())
}

/// Write TDM results as a text table to `path`.
pub fn write_tdm_paw(path: impl AsRef<Path>, entries: &[TdmPawEntry]) -> Result<()> {
    use std::io::Write;
    let mut f = std::io::BufWriter::new(std::fs::File::create(path)?);
    writeln!(f, "#{:>5} {:>5} {:>12} {:>12} {:>12} {:>12} {:>12} {:>12}",
        "iband", "jband", "E_i(eV)", "E_j(eV)", "dE(eV)", "|Tx|(D)", "|Ty|(D)", "|Tz|(D)")?;
    for e in entries {
        writeln!(f, " {:>5} {:>5} {:>12.6} {:>12.6} {:>12.6} {:>12.6} {:>12.6} {:>12.6}",
            e.iband, e.jband, e.ei, e.ej, e.de,
            e.dipole[0].norm(), e.dipole[1].norm(), e.dipole[2].norm())?;
    }
    Ok(())
}

// ---------------------------------------------------------------------------
// Helper
// ---------------------------------------------------------------------------

const AUTOA_TDM: f64 = 0.529177249;
const RYTOEV_TDM: f64 = 13.605826;

fn ngrid_from_cell(cell: &[[f64; 3]; 3], encut: f64) -> [usize; 3] {
    let anorm = [
        vec_norm(&cell[0]),
        vec_norm(&cell[1]),
        vec_norm(&cell[2]),
    ];
    std::array::from_fn(|i| {
        let cutof = ((encut / RYTOEV_TDM).sqrt()
            / (TPI / (anorm[i] / AUTOA_TDM)))
            .ceil();
        (2.0 * cutof + 1.0) as usize
    })
}

fn vec_norm(v: &[f64; 3]) -> f64 {
    (v[0] * v[0] + v[1] * v[1] + v[2] * v[2]).sqrt()
}
