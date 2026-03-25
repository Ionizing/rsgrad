// NormalCAR writer: PAW projector coefficients (beta) for all bands/k-points.

use std::fs::File;
use std::io::{BufWriter, Write};
use std::path::Path;

use anyhow::Result;
use num::complex::Complex;
use super::npy::NpzWriter;

use super::gvec::{gvectors, gvectors_gamma};
use super::nonlq::Nonlq;
use super::pawpotcar::PawPotcar;
use super::poscar::Poscar;
use super::wavecar::Wavecar;

type C128 = Complex<f64>;

const AUTOA: f64 = 0.529177249;
const RYTOEV: f64 = 13.605826;
const TPI: f64 = std::f64::consts::TAU;

/// Write a NormalCAR file from POSCAR + POTCAR + WAVECAR.
pub fn write_normalcar(
    path: impl AsRef<Path>,
    poscar: &Poscar,
    pawpot: &PawPotcar,
    wavecar: &mut Wavecar,
) -> Result<()> {
    let file = File::create(path.as_ref())?;
    let mut w = BufWriter::new(file);

    let ntypes = pawpot.len();
    let element_idx = poscar.element_idx();
    let nspin = wavecar.header.nspin;
    let nkpts = wavecar.header.nkpts;
    let nbands = wavecar.header.nbands;
    let encut = wavecar.header.encut;
    let bcell = wavecar.bcell();
    let pscut = encut;

    let lgam = detect_lgam(wavecar, poscar, pscut, &bcell);

    let lmmax_per_type: Vec<usize> = (0..ntypes)
        .map(|it| pawpot.get(it).unwrap().lmmax())
        .collect();
    let max_type = lmmax_per_type
        .iter()
        .enumerate()
        .max_by_key(|&(_, &v)| v)
        .map(|(i, _)| i)
        .unwrap_or(0);
    let lmdim = lmmax_per_type[max_type];

    let npro_tot: usize = element_idx.iter().map(|&it| lmmax_per_type[it]).sum();

    let nions_cqij: i32 = 1;
    let nrspinors: i32 = 1;
    let npro: i32 = lmdim as i32;
    let nprod: i32 = (lmdim as i32 + 3) / 4 * 4;
    let ntyp: i32 = 1;
    let natoms_max: i32 = poscar.counts[max_type] as i32;

    let qij = pawpot.get(max_type).unwrap().get_qij();
    let mut cqij: Vec<f64> = vec![0.0; lmdim * lmdim];
    for i in 0..lmdim {
        for j in 0..lmdim {
            cqij[i + j * lmdim] = qij[i][j];
        }
    }

    write_rec_i32(&mut w, &[lmdim as i32, nions_cqij, nrspinors])?;
    write_rec_f64(&mut w, &cqij)?;
    write_rec_i32(&mut w, &[nprod, npro, ntyp])?;
    write_rec_i32(&mut w, &[lmdim as i32, natoms_max])?;

    for ispin in 1..=nspin {
        for ikpt in 1..=nkpts {
            let kvec = wavecar.kpoints[ikpt - 1].kvec;
            let nonlq = Nonlq::new(poscar, pawpot, encut, &kvec, &bcell, lgam);

            for iband in 1..=nbands {
                let coeffs = wavecar.read_band_coeff(ispin, ikpt, iband)?;
                let beta = nonlq.proj(&coeffs);

                let mut rec: Vec<f64> = Vec::with_capacity(npro_tot * 2);
                for b in beta.iter().take(npro_tot) {
                    rec.push(b.re);
                    rec.push(b.im);
                }
                write_rec_f64(&mut w, &rec)?;
            }
        }
    }

    w.flush()?;
    Ok(())
}

/// Write all PAW projector coefficients (beta) to a `.npz` file.
pub fn write_cproj_npz(
    path: impl AsRef<Path>,
    poscar: &Poscar,
    pawpot: &PawPotcar,
    wavecar: &mut Wavecar,
) -> Result<()> {
    let ntypes = pawpot.len();
    let element_idx = poscar.element_idx();
    let nspin  = wavecar.header.nspin;
    let nkpts  = wavecar.header.nkpts;
    let nbands = wavecar.header.nbands;
    let encut  = wavecar.header.encut;
    let bcell  = wavecar.bcell();

    let lgam = detect_lgam(wavecar, poscar, encut, &bcell);

    let lmmax_per_type: Vec<usize> = (0..ntypes)
        .map(|it| pawpot.get(it).unwrap().lmmax())
        .collect();
    let npro_tot: usize = element_idx.iter().map(|&it| lmmax_per_type[it]).sum();

    let mut all_beta: Vec<C128> = Vec::with_capacity(nspin * nkpts * nbands * npro_tot);

    for ispin in 1..=nspin {
        for ikpt in 1..=nkpts {
            let kvec = wavecar.kpoints[ikpt - 1].kvec;
            let nonlq = Nonlq::new(poscar, pawpot, encut, &kvec, &bcell, lgam);
            for iband in 1..=nbands {
                let coeffs = wavecar.read_band_coeff(ispin, ikpt, iband)?;
                let beta = nonlq.proj(&coeffs);
                all_beta.extend_from_slice(&beta[..npro_tot]);
            }
        }
    }

    let shape = [nspin as u64, nkpts as u64, nbands as u64, npro_tot as u64];
    let mut npz = NpzWriter::create(path.as_ref())?;
    npz.write_c128("cproj", &all_beta, &shape)?;
    npz.finish()?;
    Ok(())
}

// ---------------------------------------------------------------------------
// Gamma-only detection
// ---------------------------------------------------------------------------

fn detect_lgam(wavecar: &Wavecar, poscar: &Poscar, encut: f64, bcell: &[[f64; 3]; 3]) -> bool {
    let ngrid = ngrid_from_cell(&poscar.cell, encut);
    let kvec = wavecar.kpoints[0].kvec;
    let nplw_stored = wavecar.kpoints[0].nplw;
    let nplw_full = gvectors(bcell, encut, &kvec, &ngrid).0.len();
    let nplw_gamma = gvectors_gamma(bcell, encut, &kvec, &ngrid).0.len();
    nplw_stored == nplw_gamma && nplw_stored != nplw_full
}

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

// ---------------------------------------------------------------------------
// Fortran unformatted binary helpers
// ---------------------------------------------------------------------------

fn write_rec_i32(w: &mut impl Write, data: &[i32]) -> Result<()> {
    let nbytes = (data.len() * 4) as i32;
    w.write_all(&nbytes.to_le_bytes())?;
    for &v in data {
        w.write_all(&v.to_le_bytes())?;
    }
    w.write_all(&nbytes.to_le_bytes())?;
    Ok(())
}

fn write_rec_f64(w: &mut impl Write, data: &[f64]) -> Result<()> {
    let nbytes = (data.len() * 8) as i32;
    w.write_all(&nbytes.to_le_bytes())?;
    for &v in data {
        w.write_all(&v.to_le_bytes())?;
    }
    w.write_all(&nbytes.to_le_bytes())?;
    Ok(())
}
