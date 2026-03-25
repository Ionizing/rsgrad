// AE wavefunction overlaps between two WAVECARs.

use std::collections::hash_map::Entry;

use super::npy::NpzWriter;

use anyhow::{bail, Result};
use num::complex::Complex;

use super::gvec::{gvectors, gvectors_gamma};
use super::nonlq::Nonlq;
use super::pawpotcar::PawPotcar;
use super::poscar::Poscar;
use super::wavecar::Wavecar;

type C128 = Complex<f64>;

const TPI: f64 = std::f64::consts::TAU;
const AUTOA: f64 = 0.529177249;
const RYTOEV: f64 = 13.605826;

/// One overlap entry.
#[derive(Debug)]
pub struct OverlapEntry {
    pub ispin: usize,
    pub ikpt: usize,
    pub iband: usize,
    pub jband: usize,
    /// <Φ_i^a(k) | Φ_j^b(k)>
    pub overlap: C128,
}

/// Compute overlaps <Φ_i^a | Φ_j^b> for all (ispin, ikpt, iband, jband) combinations.
#[allow(clippy::too_many_arguments)]
pub fn compute_overlap(
    poscar1: &Poscar,
    poscar2: &Poscar,
    pawpot: &PawPotcar,
    wavecar1: &mut Wavecar,
    wavecar2: &mut Wavecar,
    ispins: &[usize],
    ikpts: &[usize],
    ibands: &[usize],
    jbands: &[usize],
    lpseudo: bool,
) -> Result<Vec<OverlapEntry>> {
    let encut1 = wavecar1.header.encut;
    let encut2 = wavecar2.header.encut;
    if (encut1 - encut2).abs() > 1.0 {
        bail!("ENCUT mismatch between WAVECARs: {encut1:.1} vs {encut2:.1} eV");
    }
    if poscar1.natoms != poscar2.natoms {
        bail!("Atom count mismatch: {} vs {}", poscar1.natoms, poscar2.natoms);
    }
    if poscar1.symbols != poscar2.symbols {
        bail!("Atom types mismatch: {:?} vs {:?}", poscar1.symbols, poscar2.symbols);
    }
    if poscar1.counts != poscar2.counts {
        bail!("Atom counts mismatch: {:?} vs {:?}", poscar1.counts, poscar2.counts);
    }

    let encut = encut1;
    let bcell = wavecar1.bcell();

    let ntypes = pawpot.len();
    let qij_per_type: Vec<Vec<Vec<f64>>> = (0..ntypes)
        .map(|it| pawpot.get(it).unwrap().get_qij())
        .collect();
    let lmmax_per_type: Vec<usize> = (0..ntypes)
        .map(|it| pawpot.get(it).unwrap().lmmax())
        .collect();
    let element_idx = poscar1.element_idx();

    let mut entries = Vec::new();

    for &ispin in ispins {
        for &ikpt in ikpts {
            let kvec = wavecar1.kpoints[ikpt - 1].kvec;

            let ngrid = ngrid_from_cell(&poscar1.cell, encut);
            let nplw_stored = wavecar1.kpoints[ikpt - 1].nplw;
            let nplw_full   = gvectors(&bcell, encut, &kvec, &ngrid).0.len();
            let nplw_gamma  = gvectors_gamma(&bcell, encut, &kvec, &ngrid).0.len();
            let lgam = nplw_stored == nplw_gamma && nplw_stored != nplw_full;

            let nonlq1 = Nonlq::new(poscar1, pawpot, encut, &kvec, &bcell, lgam);
            let nonlq2 = Nonlq::new(poscar2, pawpot, encut, &kvec, &bcell, lgam);
            let nplw = nonlq1.nplw;

            let mut cache1: std::collections::HashMap<usize, Vec<C128>> = Default::default();
            for &b in ibands {
                if let Entry::Vacant(e) = cache1.entry(b) {
                    let c = wavecar1.read_band_coeff(ispin, ikpt, b)?;
                    e.insert(c);
                }
            }
            let mut cache2: std::collections::HashMap<usize, Vec<C128>> = Default::default();
            for &b in jbands {
                if let Entry::Vacant(e) = cache2.entry(b) {
                    let c = wavecar2.read_band_coeff(ispin, ikpt, b)?;
                    e.insert(c);
                }
            }

            let mut beta2_cache: std::collections::HashMap<usize, Vec<C128>> = Default::default();
            if !lpseudo {
                for &jb in jbands {
                    if let Entry::Vacant(e) = beta2_cache.entry(jb) {
                        let bj = nonlq2.proj(&cache2[&jb]);
                        e.insert(bj);
                    }
                }
            }

            for &ib in ibands {
                let ci = &cache1[&ib];

                let beta1_i = if !lpseudo {
                    Some(nonlq1.proj(ci))
                } else {
                    None
                };

                for &jb in jbands {
                    let cj = &cache2[&jb];

                    let ps_ovlp: C128 = {
                        let raw: C128 = (0..nplw).map(|ig| ci[ig].conj() * cj[ig]).sum();
                        if lgam { C128::new(raw.re, 0.0) } else { raw }
                    };

                    let paw_corr = if !lpseudo {
                        let bi = beta1_i.as_ref().unwrap();
                        let bj = &beta2_cache[&jb];

                        let mut corr = C128::new(0.0, 0.0);
                        let mut nproj = 0usize;
                        for &ntype in &element_idx {
                            let lmmax = lmmax_per_type[ntype];
                            let qij = &qij_per_type[ntype];
                            for m in 0..lmmax {
                                let mut tmp = C128::new(0.0, 0.0);
                                for n in 0..lmmax {
                                    tmp += qij[m][n] * bj[nproj + n];
                                }
                                corr += bi[nproj + m].conj() * tmp;
                            }
                            nproj += lmmax;
                        }
                        if lgam { C128::new(corr.re, 0.0) } else { corr }
                    } else {
                        C128::new(0.0, 0.0)
                    };

                    entries.push(OverlapEntry {
                        ispin,
                        ikpt,
                        iband: ib,
                        jband: jb,
                        overlap: ps_ovlp + paw_corr,
                    });
                }
            }
        }
    }

    Ok(entries)
}

/// Write overlap results to a `.npz` file.
pub fn write_overlap_npz(
    path: impl AsRef<std::path::Path>,
    entries: &[OverlapEntry],
) -> Result<()> {
    let n = entries.len();
    let ispin:   Vec<i64>  = entries.iter().map(|e| e.ispin as i64).collect();
    let ikpt:    Vec<i64>  = entries.iter().map(|e| e.ikpt  as i64).collect();
    let iband:   Vec<i64>  = entries.iter().map(|e| e.iband as i64).collect();
    let jband:   Vec<i64>  = entries.iter().map(|e| e.jband as i64).collect();
    let overlap: Vec<C128> = entries.iter().map(|e| e.overlap).collect();

    let mut npz = NpzWriter::create(path.as_ref())?;
    npz.write_i64("ispin",   &ispin,   &[n as u64])?;
    npz.write_i64("ikpt",    &ikpt,    &[n as u64])?;
    npz.write_i64("iband",   &iband,   &[n as u64])?;
    npz.write_i64("jband",   &jband,   &[n as u64])?;
    npz.write_c128("overlap", &overlap, &[n as u64])?;
    npz.finish()?;
    Ok(())
}

/// Write overlap results as a human-readable text table.
pub fn write_overlap(
    path: impl AsRef<std::path::Path>,
    entries: &[OverlapEntry],
) -> Result<()> {
    use std::io::Write;
    let mut f = std::io::BufWriter::new(std::fs::File::create(path)?);
    writeln!(
        f,
        "#{:>4} {:>4} {:>5} {:>5} {:>18} {:>18} {:>14}",
        "spin", "kpt", "iband", "jband", "Re(<i|j>)", "Im(<i|j>)", "|<i|j>|"
    )?;
    for e in entries {
        writeln!(
            f,
            " {:>4} {:>4} {:>5} {:>5} {:>18.10} {:>18.10} {:>14.10}",
            e.ispin, e.ikpt, e.iband, e.jband,
            e.overlap.re, e.overlap.im, e.overlap.norm()
        )?;
    }
    Ok(())
}

// ---------------------------------------------------------------------------
// Helpers
// ---------------------------------------------------------------------------

fn ngrid_from_cell(cell: &[[f64; 3]; 3], encut: f64) -> [usize; 3] {
    let anorm = [
        vec_norm(&cell[0]),
        vec_norm(&cell[1]),
        vec_norm(&cell[2]),
    ];
    std::array::from_fn(|i| {
        let cutof = ((encut / RYTOEV).sqrt() / (TPI / (anorm[i] / AUTOA))).ceil();
        (2.0 * cutof + 1.0) as usize
    })
}

fn vec_norm(v: &[f64; 3]) -> f64 {
    (v[0] * v[0] + v[1] * v[1] + v[2] * v[2]).sqrt()
}
