// G-vector generator: enumerate G-vectors satisfying (G+k)²/2 < ENCUT.
//
// This follows the exact ordering used by VASP/VaspBandUnfolding,
// which uses np.meshgrid(fz, fy, fx, indexing='ij') and then reshapes to
// (gx, gy, gz) — making gx fastest-varying, gz slowest.

use super::wavecar::{matvec_row, Wavecar};

// Physical constants
const AUTOA: f64 = 0.529177249;
const RYTOEV: f64 = 13.605826;
const HSQDTM: f64 = RYTOEV * AUTOA * AUTOA;
const TPI: f64 = std::f64::consts::TAU;

/// Generate G-vectors (as integer triplets [gx, gy, gz]) and their
/// Cartesian counterparts `(G+k)` for a given k-point.
pub fn gvectors(
    bcell: &[[f64; 3]; 3],
    encut: f64,
    kvec: &[f64; 3],
    ngrid: &[usize; 3],
) -> (Vec<[i32; 3]>, Vec<[f64; 3]>, Vec<f64>) {
    let [nx, ny, nz] = *ngrid;

    let fx: Vec<i32> = fft_indices(nx);
    let fy: Vec<i32> = fft_indices(ny);
    let fz: Vec<i32> = fft_indices(nz);

    let mut gvec = Vec::new();
    let mut gk_cart = Vec::new();
    let mut glen = Vec::new();

    for &gz in &fz {
        for &gy in &fy {
            for &gx in &fx {
                let kfrac = [
                    gx as f64 + kvec[0],
                    gy as f64 + kvec[1],
                    gz as f64 + kvec[2],
                ];
                let gk = matvec_row(bcell, &[kfrac[0] * TPI, kfrac[1] * TPI, kfrac[2] * TPI]);
                let g2 = gk[0] * gk[0] + gk[1] * gk[1] + gk[2] * gk[2];
                let ke = HSQDTM * g2;
                if ke < encut {
                    gvec.push([gx, gy, gz]);
                    gk_cart.push(gk);
                    glen.push(g2.sqrt());
                }
            }
        }
    }

    (gvec, gk_cart, glen)
}

/// FFT-ordered integer indices for dimension n:
/// [0, 1, ..., n/2, -(n/2 - 1), ..., -1]
fn fft_indices(n: usize) -> Vec<i32> {
    let mut v: Vec<i32> = (0..n as i32).collect();
    let half = n / 2 + 1;
    for x in v[half..].iter_mut() {
        *x -= n as i32;
    }
    v
}

/// Generate gamma-only G-vectors (gamma_half='x' convention).
pub fn gvectors_gamma(
    bcell: &[[f64; 3]; 3],
    encut: f64,
    kvec: &[f64; 3],
    ngrid: &[usize; 3],
) -> (Vec<[i32; 3]>, Vec<[f64; 3]>, Vec<f64>) {
    let [nx, ny, nz] = *ngrid;

    let fx: Vec<i32> = (0..=(nx / 2) as i32).collect();
    let fy: Vec<i32> = fft_indices(ny);
    let fz: Vec<i32> = fft_indices(nz);

    let mut gvec = Vec::new();
    let mut gk_cart = Vec::new();
    let mut glen = Vec::new();

    for &gz in &fz {
        for &gy in &fy {
            for &gx in &fx {
                if !(gx > 0 || (gx == 0 && gy > 0) || (gx == 0 && gy == 0 && gz >= 0)) {
                    continue;
                }
                let kfrac = [
                    gx as f64 + kvec[0],
                    gy as f64 + kvec[1],
                    gz as f64 + kvec[2],
                ];
                let gk = matvec_row(bcell, &[kfrac[0] * TPI, kfrac[1] * TPI, kfrac[2] * TPI]);
                let g2 = gk[0] * gk[0] + gk[1] * gk[1] + gk[2] * gk[2];
                let ke = HSQDTM * g2;
                if ke < encut {
                    gvec.push([gx, gy, gz]);
                    gk_cart.push(gk);
                    glen.push(g2.sqrt());
                }
            }
        }
    }

    (gvec, gk_cart, glen)
}

/// Use the Wavecar's bcell and ngrid to generate G-vectors for a k-point.
pub fn gvectors_from_wavecar(
    wfc: &Wavecar,
    ikpt: usize, // 1-indexed
) -> (Vec<[i32; 3]>, Vec<[f64; 3]>, Vec<f64>) {
    let bcell = wfc.bcell();
    let kvec = wfc.kpoints[ikpt - 1].kvec;
    gvectors(&bcell, wfc.header.encut, &kvec, &wfc.header.ngrid)
}

#[cfg(test)]
mod tests {
    use super::*;

    #[test]
    fn test_gvectors_count_matches_wavecar() {
        let wfc = Wavecar::from_file("tests/pawpot/projectors_lreal_false/WAVECAR").unwrap();
        for ikpt in 1..=3 {
            let (gvec, _, _) = gvectors_from_wavecar(&wfc, ikpt);
            let expected = wfc.kpoints[ikpt - 1].nplw;
            assert_eq!(
                gvec.len(),
                expected,
                "ikpt={ikpt}: got {} G-vectors, expected {expected}",
                gvec.len()
            );
        }
    }
}
