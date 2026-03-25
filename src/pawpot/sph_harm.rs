// Real spherical harmonics following the VaspBandUnfolding / scipy convention.
//
// Convention for real Y_l^m (ordering: m = -l, ..., 0, ..., +l):
//   m < 0: Y_R = (-1)^|m| * sqrt(2) * Im(Y_C^|m|)
//   m = 0: Y_R = Re(Y_C^0)
//   m > 0: Y_R = (-1)^m  * sqrt(2) * Re(Y_C^m)
//
// where Y_C^m is the complex SH with Condon-Shortley phase.

use std::f64::consts::{PI, SQRT_2};
use num::complex::Complex64;

const EPSILON: f64 = 1e-10;

/// Convert Cartesian coordinates to spherical (r, phi, theta).
/// phi: azimuthal angle in [0, 2π], theta: polar angle in [0, π].
pub fn cart2sph(xyz: &[[f64; 3]]) -> (Vec<f64>, Vec<f64>, Vec<f64>) {
    let n = xyz.len();
    let mut r_out = vec![0.0; n];
    let mut phi_out = vec![0.0; n];
    let mut theta_out = vec![0.0; n];

    for i in 0..n {
        let [x, y, z] = xyz[i];
        let ri = (x * x + y * y + z * z).sqrt().max(EPSILON);
        let mut p = y.atan2(x);
        if p < 0.0 {
            p += 2.0 * PI;
        }
        r_out[i] = ri;
        phi_out[i] = p;
        theta_out[i] = (z / ri).clamp(-1.0, 1.0).acos();
    }
    (r_out, phi_out, theta_out)
}

/// Associated Legendre polynomial P_l^m(x) with Condon-Shortley phase, m >= 0.
fn assoc_legendre(l: usize, m: usize, x: f64) -> f64 {
    debug_assert!(m <= l);

    let mut pmm = 1.0_f64;
    let somx2 = ((1.0 - x) * (1.0 + x)).sqrt();
    let mut fact = 1.0_f64;
    for _ in 0..m {
        pmm *= -fact * somx2;
        fact += 2.0;
    }
    if l == m {
        return pmm;
    }

    let mut pmmp1 = x * (2 * m + 1) as f64 * pmm;
    if l == m + 1 {
        return pmmp1;
    }

    let mut pll = 0.0_f64;
    for ll in (m + 2)..=l {
        pll = (x * (2 * ll - 1) as f64 * pmmp1 - (ll + m - 1) as f64 * pmm) / (ll - m) as f64;
        pmm = pmmp1;
        pmmp1 = pll;
    }
    pll
}

fn factorial(n: usize) -> f64 {
    (1..=n).fold(1.0_f64, |acc, x| acc * x as f64)
}

/// Complex spherical harmonic Y_l^m for m >= 0 (includes Condon-Shortley phase).
fn ylm_positive_m(l: usize, m: usize, theta: f64, phi: f64) -> Complex64 {
    let plm = assoc_legendre(l, m, theta.cos());
    let norm = ((2 * l + 1) as f64 / (4.0 * PI) * factorial(l - m) / factorial(l + m)).sqrt();
    let (sin_mphi, cos_mphi) = (m as f64 * phi).sin_cos();
    norm * plm * Complex64::new(cos_mphi, sin_mphi)
}

/// Real spherical harmonics for all m = -l..=l at multiple points.
///
/// Returns a `Vec` of length `n_points`, each inner `Vec` has length `2*l+1`.
/// Ordering: m = -l, ..., 0, ..., l.
pub fn sph_r(xyz: &[[f64; 3]], l: usize) -> Vec<Vec<f64>> {
    let (_, phi_arr, theta_arr) = cart2sph(xyz);
    let n = xyz.len();
    let tlp1 = 2 * l + 1;

    let mut result = vec![vec![0.0; tlp1]; n];

    for i in 0..n {
        let theta = theta_arr[i];
        let phi = phi_arr[i];

        // Complex Y_l^m for m = 0..=l
        let ylm_c: Vec<Complex64> = (0..=l)
            .map(|m| ylm_positive_m(l, m, theta, phi))
            .collect();

        for (ii, m_signed) in (-(l as i64)..=(l as i64)).enumerate() {
            result[i][ii] = if m_signed < 0 {
                let am = (-m_signed) as usize;
                let factor = if am % 2 == 0 { 1.0 } else { -1.0 };
                factor * SQRT_2 * ylm_c[am].im
            } else if m_signed == 0 {
                ylm_c[0].re
            } else {
                let m = m_signed as usize;
                let factor = if m % 2 == 0 { 1.0 } else { -1.0 };
                factor * SQRT_2 * ylm_c[m].re
            };
        }
    }
    result
}

#[cfg(test)]
mod tests {
    use super::*;

    #[test]
    fn test_sph_r_l0_constant() {
        let xyz = vec![
            [1.0, 0.0, 0.0],
            [0.0, 1.0, 0.0],
            [0.0, 0.0, 1.0],
            [1.0, 1.0, 1.0],
        ];
        let expected = 1.0 / (4.0 * PI).sqrt(); // 0.28209479
        let result = sph_r(&xyz, 0);
        for row in &result {
            assert!((row[0] - expected).abs() < 1e-10);
        }
    }

    #[test]
    fn test_sph_r_l1_axes() {
        let tol = 1e-9;
        // At (1,0,0): m=-1 -> 0, m=0 -> 0, m=1 -> sqrt(3/(4π))
        let result_x = sph_r(&[[1.0, 0.0, 0.0]], 1);
        let y10 = (3.0 / (4.0 * PI)).sqrt(); // ≈ 0.48860251
        assert!(result_x[0][0].abs() < tol); // m=-1
        assert!(result_x[0][1].abs() < tol); // m=0
        assert!((result_x[0][2] - y10).abs() < tol); // m=1

        // At (0,1,0): m=-1 -> sqrt(3/(4π)), m=0 -> 0, m=1 -> 0
        let result_y = sph_r(&[[0.0, 1.0, 0.0]], 1);
        assert!((result_y[0][0] - y10).abs() < tol); // m=-1
        assert!(result_y[0][1].abs() < tol); // m=0
        assert!(result_y[0][2].abs() < tol); // m=1

        // At (0,0,1): m=-1 -> 0, m=0 -> sqrt(3/(4π)), m=1 -> 0
        let result_z = sph_r(&[[0.0, 0.0, 1.0]], 1);
        assert!(result_z[0][0].abs() < tol); // m=-1
        assert!((result_z[0][1] - y10).abs() < tol); // m=0
        assert!(result_z[0][2].abs() < tol); // m=1
    }

    #[test]
    fn test_sph_r_normalization() {
        // Monte Carlo: average of Y^2 over uniform sphere ≈ 1/(4π)
        let n = 5000_usize;
        let mut state = 0x123456789abcdef0_u64;
        let mut rng = || -> f64 {
            state ^= state << 13;
            state ^= state >> 7;
            state ^= state << 17;
            (state >> 11) as f64 / (1u64 << 53) as f64
        };
        let xyz: Vec<[f64; 3]> = (0..n)
            .map(|_| {
                let cos_t = rng() * 2.0 - 1.0;
                let phi = rng() * 2.0 * PI;
                let sin_t = (1.0 - cos_t * cos_t).sqrt();
                [sin_t * phi.cos(), sin_t * phi.sin(), cos_t]
            })
            .collect();

        for l in 0usize..=2 {
            let vals = sph_r(&xyz, l);
            for m_idx in 0..(2 * l + 1) {
                let sum_sq: f64 = vals.iter().map(|row| row[m_idx] * row[m_idx]).sum();
                // integral Y^2 dΩ = 1, so average = 1/(4π)
                let expected = n as f64 / (4.0 * PI);
                let rel_err = (sum_sq - expected).abs() / expected;
                assert!(
                    rel_err < 0.05,
                    "l={}, m_idx={}: sum_sq={:.4}, expected={:.4}, rel_err={:.4}",
                    l, m_idx, sum_sq, expected, rel_err
                );
            }
        }
    }
}
