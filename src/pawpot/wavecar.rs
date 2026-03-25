// WAVECAR reader for VASP plane-wave wavefunction files.
//
// Binary format (Fortran unformatted, record-based):
//   Record 0: recl (f64), nspin (f64), rtag (f64)
//   Record 1: nkpts, nbands, encut, cell(9), efermi  [13 f64]
//   For each spin/kpt: k-point header + band coefficient records.
//
// rtag = 45200 → complex64 coefficients
// rtag = 45210 → complex128 coefficients

use std::convert::TryInto;

use anyhow::{bail, Context, Result};
use num::complex::Complex;
use std::fs::File;
use std::io::{Read, Seek, SeekFrom};
use std::path::Path;

pub type Complex64 = Complex<f32>;
pub type Complex128 = Complex<f64>;

/// Wavefunction precision tag.
#[derive(Debug, Clone, Copy, PartialEq)]
pub enum WfPrec {
    Single,  // rtag = 45200, complex64
    Double,  // rtag = 45210, complex128
}

/// WAVECAR header information.
#[derive(Debug, Clone)]
pub struct WavecarHeader {
    pub recl: u64,
    pub nspin: usize,
    pub rtag: usize,
    pub prec: WfPrec,
    pub nkpts: usize,
    pub nbands: usize,
    pub encut: f64,
    /// Real-space lattice vectors, cell[i] = [x, y, z] in Å.
    pub cell: [[f64; 3]; 3],
    pub efermi: f64,
    /// Minimum FFT grid size [nx, ny, nz].
    pub ngrid: [usize; 3],
}

/// Per-k-point data (k-vector and per-band metadata).
#[derive(Debug, Clone)]
pub struct KpointInfo {
    pub nplw: usize,
    pub kvec: [f64; 3],
    /// Eigenvalues in eV.
    pub eigenvalues: Vec<f64>,
    /// Occupations.
    pub occupations: Vec<f64>,
}

/// WAVECAR reader.
pub struct Wavecar {
    pub header: WavecarHeader,
    pub kpoints: Vec<KpointInfo>,
    file: File,
}

// Physical constants (from VASP)
const AUTOA: f64 = 0.529177249;
const RYTOEV: f64 = 13.605826;
const HSQDTM: f64 = RYTOEV * AUTOA * AUTOA;
const TPI: f64 = std::f64::consts::TAU;

impl Wavecar {
    pub fn from_file(path: impl AsRef<Path>) -> Result<Self> {
        let mut file = File::open(path.as_ref())
            .with_context(|| format!("opening {}", path.as_ref().display()))?;

        // Record 0: recl, nspin, rtag
        let mut buf3 = [0u8; 24];
        file.read_exact(&mut buf3)?;
        let recl = f64::from_le_bytes(buf3[0..8].try_into().unwrap()) as u64;
        let nspin = f64::from_le_bytes(buf3[8..16].try_into().unwrap()) as usize;
        let rtag = f64::from_le_bytes(buf3[16..24].try_into().unwrap()) as usize;

        let prec = match rtag {
            45200 => WfPrec::Single,
            45210 => WfPrec::Double,
            _ => bail!("unsupported WAVECAR rtag: {rtag}"),
        };

        // Record 1: nkpts, nbands, encut, cell(9), efermi
        file.seek(SeekFrom::Start(recl))?;
        let hdr = read_f64s(&mut file, 13)?;
        let nkpts = hdr[0] as usize;
        let nbands = hdr[1] as usize;
        let encut = hdr[2];
        let mut cell = [[0.0f64; 3]; 3];
        for i in 0..3 {
            cell[i] = [hdr[3 + 3 * i], hdr[3 + 3 * i + 1], hdr[3 + 3 * i + 2]];
        }
        let efermi = hdr[12];

        // Compute minimum FFT grid
        let ngrid = compute_ngrid(&cell, encut);

        let header = WavecarHeader {
            recl,
            nspin,
            rtag,
            prec,
            nkpts,
            nbands,
            encut,
            cell,
            efermi,
            ngrid,
        };

        // Read k-point info for spin 1 only (spin-polarised not needed yet)
        let mut kpoints = Vec::with_capacity(nkpts);
        for ikpt in 0..nkpts {
            let rec = 2 + ikpt * (nbands + 1); // 0-indexed record number for kpt header
            file.seek(SeekFrom::Start(rec as u64 * recl))?;
            let dump = read_f64s(&mut file, 4 + 3 * nbands)?;
            let nplw = dump[0] as usize;
            let kvec = [dump[1], dump[2], dump[3]];
            let mut eigenvalues = Vec::with_capacity(nbands);
            let mut occupations = Vec::with_capacity(nbands);
            for ib in 0..nbands {
                eigenvalues.push(dump[4 + 3 * ib]);
                occupations.push(dump[4 + 3 * ib + 2]);
            }
            kpoints.push(KpointInfo { nplw, kvec, eigenvalues, occupations });
        }

        Ok(Wavecar { header, kpoints, file })
    }

    /// Record number (0-indexed) for ispin (1-indexed), ikpt (1-indexed), iband (1-indexed).
    fn where_rec(&self, ispin: usize, ikpt: usize, iband: usize) -> u64 {
        (2 + (ispin - 1) * self.header.nkpts * (self.header.nbands + 1)
            + (ikpt - 1) * (self.header.nbands + 1)
            + iband) as u64
    }

    /// Read plane-wave coefficients as complex128 for the given state.
    /// ispin, ikpt, iband are 1-indexed.
    pub fn read_band_coeff(&mut self, ispin: usize, ikpt: usize, iband: usize) -> Result<Vec<Complex128>> {
        let rec = self.where_rec(ispin, ikpt, iband);
        self.file.seek(SeekFrom::Start(rec * self.header.recl))?;
        let nplw = self.kpoints[ikpt - 1].nplw;

        let coeffs = match self.header.prec {
            WfPrec::Single => {
                let mut buf = vec![0u8; nplw * 8];
                self.file.read_exact(&mut buf)?;
                buf.chunks_exact(8)
                    .map(|c| {
                        let re = f32::from_le_bytes(c[0..4].try_into().unwrap()) as f64;
                        let im = f32::from_le_bytes(c[4..8].try_into().unwrap()) as f64;
                        Complex128::new(re, im)
                    })
                    .collect()
            }
            WfPrec::Double => {
                let mut buf = vec![0u8; nplw * 16];
                self.file.read_exact(&mut buf)?;
                buf.chunks_exact(16)
                    .map(|c| {
                        let re = f64::from_le_bytes(c[0..8].try_into().unwrap());
                        let im = f64::from_le_bytes(c[8..16].try_into().unwrap());
                        Complex128::new(re, im)
                    })
                    .collect()
            }
        };
        Ok(coeffs)
    }

    /// Reciprocal lattice vectors B = inv(A)^T, in Å^{-1}.
    pub fn bcell(&self) -> [[f64; 3]; 3] {
        let c = &self.header.cell;
        let vol = det3(c);
        let b0 = cross(&c[1], &c[2]);
        let b1 = cross(&c[2], &c[0]);
        let b2 = cross(&c[0], &c[1]);
        [
            [b0[0] / vol, b0[1] / vol, b0[2] / vol],
            [b1[0] / vol, b1[1] / vol, b1[2] / vol],
            [b2[0] / vol, b2[1] / vol, b2[2] / vol],
        ]
    }
}

/// Read n f64 values from reader.
fn read_f64s(f: &mut impl Read, n: usize) -> Result<Vec<f64>> {
    let mut buf = vec![0u8; n * 8];
    f.read_exact(&mut buf)?;
    Ok(buf
        .chunks_exact(8)
        .map(|c| f64::from_le_bytes(c.try_into().unwrap()))
        .collect())
}

/// Compute minimum FFT grid from cell and energy cutoff.
fn compute_ngrid(cell: &[[f64; 3]; 3], encut: f64) -> [usize; 3] {
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

fn cross(a: &[f64; 3], b: &[f64; 3]) -> [f64; 3] {
    [
        a[1] * b[2] - a[2] * b[1],
        a[2] * b[0] - a[0] * b[2],
        a[0] * b[1] - a[1] * b[0],
    ]
}

fn det3(m: &[[f64; 3]; 3]) -> f64 {
    let [a, b, c] = m;
    a[0] * (b[1] * c[2] - b[2] * c[1])
        - a[1] * (b[0] * c[2] - b[2] * c[0])
        + a[2] * (b[0] * c[1] - b[1] * c[0])
}

/// Compute HSQDTM * |dot(kfrac + kvec, TPI*Bcell)|² in eV.
pub fn kinetic_energy(kfrac: &[f64; 3], kvec: &[f64; 3], bcell: &[[f64; 3]; 3]) -> f64 {
    let gk = matvec_mul(bcell, &[
        (kfrac[0] + kvec[0]) * TPI,
        (kfrac[1] + kvec[1]) * TPI,
        (kfrac[2] + kvec[2]) * TPI,
    ]);
    HSQDTM * (gk[0] * gk[0] + gk[1] * gk[1] + gk[2] * gk[2])
}

/// y = M^T x  (treating rows of M as basis vectors, so Gk = (k+G) in Cartesian)
/// In numpy: np.dot(kgrid + kvec, TPI*Bcell) where Bcell is row-major.
/// We want the row-vector * matrix product, i.e., sum_j x[j] * M[j][i].
pub fn matvec_row(m: &[[f64; 3]; 3], x: &[f64; 3]) -> [f64; 3] {
    [
        x[0] * m[0][0] + x[1] * m[1][0] + x[2] * m[2][0],
        x[0] * m[0][1] + x[1] * m[1][1] + x[2] * m[2][1],
        x[0] * m[0][2] + x[1] * m[1][2] + x[2] * m[2][2],
    ]
}

/// y = M x (column-vector matrix multiply)
pub fn matvec_mul(m: &[[f64; 3]; 3], x: &[f64; 3]) -> [f64; 3] {
    [
        m[0][0] * x[0] + m[0][1] * x[1] + m[0][2] * x[2],
        m[1][0] * x[0] + m[1][1] * x[1] + m[1][2] * x[2],
        m[2][0] * x[0] + m[2][1] * x[1] + m[2][2] * x[2],
    ]
}

#[cfg(test)]
mod tests {
    use super::*;

    #[test]
    fn test_wavecar_header() {
        let wfc = Wavecar::from_file("tests/pawpot/projectors_lreal_false/WAVECAR").unwrap();
        let h = &wfc.header;
        assert_eq!(h.nspin, 1);
        assert_eq!(h.rtag, 45200);
        assert_eq!(h.nkpts, 3);
        assert_eq!(h.nbands, 16);
        assert!((h.encut - 323.36).abs() < 0.01, "encut = {}", h.encut);
        // ngrid should be [11, 11, 61]
        assert_eq!(h.ngrid, [11, 11, 61]);
    }

    #[test]
    fn test_wavecar_kpoints() {
        let wfc = Wavecar::from_file("tests/pawpot/projectors_lreal_false/WAVECAR").unwrap();
        assert_eq!(wfc.kpoints[0].nplw, 2363);
        assert_eq!(wfc.kpoints[1].nplw, 2324);
        assert_eq!(wfc.kpoints[2].nplw, 2295);

        // kvec[0] ≈ (0, 0, 0)
        assert!(wfc.kpoints[0].kvec[0].abs() < 1e-6);
        // kvec[1] ≈ (1/3, 0, 0)
        assert!((wfc.kpoints[1].kvec[0] - 1.0 / 3.0).abs() < 1e-6);
    }
}
