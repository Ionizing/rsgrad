//! **This file is adapted from [scirs2-special](https://crates.io/crates/scirs2-special)
//! under MIT license.**
//!
//! Spherical Bessel functions
//!
//! This module provides implementations of spherical Bessel functions
//! with enhanced numerical stability for both small and large arguments.
//!
//! Spherical Bessel functions are solutions to the differential equation:
//! x² d²y/dx² + 2x dy/dx + [x² - n(n+1)]y = 0
//!
//! Functions included in this module:
//! - spherical_jn(n, x): Spherical Bessel function of the first kind
//! - spherical_yn(n, x): Spherical Bessel function of the second kind
//! - spherical_jn_scaled(n, x): Scaled spherical Bessel function of the first kind (for large arguments)
//! - spherical_yn_scaled(n, x): Scaled spherical Bessel function of the second kind (for large arguments)

use num::{Float, FromPrimitive};
use std::fmt::Debug;

/// Helper function for small argument series expansion of spherical Bessel functions
fn small_arg_series_jn<F: Float + FromPrimitive + Debug>(n: i32, x: F) -> F {
    let _n_f = F::from(n).unwrap();
    let x_sq = x * x;

    let mut factorial = F::one();
    for i in 1..=n {
        factorial = factorial * F::from(2 * i + 1).unwrap();
    }

    let mut term = F::from(1.0).unwrap();
    let mut series = term;

    let terms_to_compute = if n < 5 { 4 } else { 3 };

    for i in 1..=terms_to_compute {
        let denom = F::from(2.0 * i as f64).unwrap() * F::from((2 * n + 1 + 2 * i) as f64).unwrap();
        term = term * x_sq.neg() / denom;
        series = series + term;

        if term.abs() < F::from(1e-15).unwrap() * series.abs() {
            break;
        }
    }

    let mut x_pow_n = F::one();
    for _ in 0..n {
        x_pow_n = x_pow_n * x;
    }

    x_pow_n / factorial * series
}

/// Spherical Bessel function of the first kind with enhanced stability.
pub fn spherical_jn<F: Float + FromPrimitive + Debug>(n: i32, x: F) -> F {
    if n < 0 {
        panic!("Order n must be non-negative");
    }

    if x == F::zero() {
        if n == 0 {
            return F::one();
        } else {
            return F::zero();
        }
    }

    if n == 0 {
        if x.abs() < F::from(0.01).unwrap() {
            let x2 = x * x;
            return F::one() - x2 / F::from(6.0).unwrap() + x2 * x2 / F::from(120.0).unwrap();
        } else {
            return x.sin() / x;
        }
    } else if n == 1 {
        if x.abs() < F::from(0.01).unwrap() {
            let x2 = x * x;
            return x / F::from(3.0).unwrap() - x * x2 / F::from(30.0).unwrap();
        } else {
            return (x.sin() / x - x.cos()) / x;
        }
    }

    if x.abs() < F::from(0.1).unwrap() * (F::from(n).unwrap() + F::one()) {
        return small_arg_series_jn(n, x);
    }

    if x > F::from(n).unwrap() * F::from(10.0).unwrap() {
        let scaling = x.sin();
        return spherical_jn_scaled(n, x) * scaling / x;
    }

    let max_n = n.min(1000);

    let mut j_n_minus_2 = if x.abs() < F::from(0.01).unwrap() {
        let x2 = x * x;
        F::one() - x2 / F::from(6.0).unwrap() + x2 * x2 / F::from(120.0).unwrap()
    } else {
        x.sin() / x
    };
    let mut j_n_minus_1 = if x.abs() < F::from(0.01).unwrap() {
        let x2 = x * x;
        x / F::from(3.0).unwrap() - x * x2 / F::from(30.0).unwrap()
    } else {
        (x.sin() / x - x.cos()) / x
    };

    for k in 2..=max_n {
        let j_n = F::from(2.0 * k as f64 - 1.0).unwrap() / x * j_n_minus_1 - j_n_minus_2;
        j_n_minus_2 = j_n_minus_1;
        j_n_minus_1 = j_n;
    }

    j_n_minus_1
}

/// Spherical Bessel function of the second kind with enhanced stability.
pub fn spherical_yn<F: Float + FromPrimitive + Debug>(n: i32, x: F) -> F {
    if n < 0 {
        panic!("Order n must be non-negative");
    }

    if x <= F::zero() {
        return F::neg_infinity();
    }

    if n == 0 {
        return -x.cos() / x;
    } else if n == 1 {
        return -(x.cos() / x + x.sin()) / x;
    }

    if x > F::from(n).unwrap() * F::from(10.0).unwrap() {
        let scaling = -x.cos();
        return spherical_yn_scaled(n, x) * scaling / x;
    }

    let max_n = n.min(1000);

    let mut y_n_minus_2 = -x.cos() / x;
    let mut y_n_minus_1 = -(x.cos() / x + x.sin()) / x;

    for k in 2..=max_n {
        let y_n = F::from(2.0 * k as f64 - 1.0).unwrap() / x * y_n_minus_1 - y_n_minus_2;
        y_n_minus_2 = y_n_minus_1;
        y_n_minus_1 = y_n;
    }

    y_n_minus_1
}

/// Scaled spherical Bessel function of the first kind.
pub fn spherical_jn_scaled<F: Float + FromPrimitive + Debug>(n: i32, x: F) -> F {
    if n < 0 {
        panic!("Order n must be non-negative");
    }

    if n == 0 {
        if x == F::zero() {
            return F::one();
        }

        if x > F::from(10.0).unwrap() {
            let x_sq = x * x;
            return F::one() - F::one() / (F::from(2.0).unwrap() * x_sq);
        } else {
            return F::one();
        }
    }

    if n == 1 {
        if x == F::zero() {
            return F::zero();
        }

        if x > F::from(10.0).unwrap() {
            let x_sq = x * x;
            return F::one() - F::from(3.0).unwrap() / (F::from(2.0).unwrap() * x_sq);
        } else {
            let x2 = x * x;
            return F::one() - F::from(2.0).unwrap() / F::from(3.0).unwrap() * x2;
        }
    }

    if x < F::from(5.0).unwrap() {
        if x == F::zero() {
            return F::zero();
        }

        let j_n = spherical_jn(n, x);
        return j_n * x / x.sin();
    }

    if x > F::from(n * n).unwrap() || x > F::from(1000.0).unwrap() {
        let x_sq = x * x;
        let n_f = F::from(n).unwrap();

        let mut factor = F::one();

        factor = factor - n_f * (n_f + F::one()) / (F::from(2.0).unwrap() * x_sq);

        if x > F::from(100.0).unwrap() {
            let term2 = n_f * (n_f + F::one()) * (n_f * (n_f + F::one()) - F::from(2.0).unwrap())
                / (F::from(8.0).unwrap() * x_sq * x_sq);
            factor = factor + term2;
        }

        return factor;
    }

    let safe_n = n.min(50);

    let n_max = (n * 2).min(100);

    let mut j_n_plus_1 = F::from(1e-100).unwrap();
    let mut j_n = F::from(1e-100).unwrap();

    for k in (0..=n_max).rev() {
        let j_n_minus_1 = F::from(2.0 * k as f64 + 1.0).unwrap() / x * j_n - j_n_plus_1;
        j_n_plus_1 = j_n;
        j_n = j_n_minus_1;

        if j_n.abs() > F::from(1e50).unwrap() {
            let scale = F::from(1e-50).unwrap();
            j_n = j_n * scale;
            j_n_plus_1 = j_n_plus_1 * scale;
        }
    }

    let j_0_scaled = F::one();

    let scale = j_0_scaled / j_n;

    j_n = j_0_scaled;
    j_n_plus_1 = j_n_plus_1 * scale;

    for k in 0..safe_n {
        let j_n_plus_1_new = F::from(2.0 * k as f64 + 1.0).unwrap() / x * j_n - j_n_plus_1;
        j_n_plus_1 = j_n;
        j_n = j_n_plus_1_new;
    }

    j_n
}

/// Scaled spherical Bessel function of the second kind.
pub fn spherical_yn_scaled<F: Float + FromPrimitive + Debug>(n: i32, x: F) -> F {
    if n < 0 {
        panic!("Order n must be non-negative");
    }

    if x <= F::zero() {
        return F::neg_infinity();
    }

    if n == 0 {
        if x > F::from(10.0).unwrap() {
            let x_sq = x * x;
            return F::one() - F::one() / (F::from(2.0).unwrap() * x_sq);
        } else {
            return -x.cos() / x * x / (-x.cos());
        }
    }

    if n == 1 {
        if x > F::from(10.0).unwrap() {
            let x_sq = x * x;
            return -F::one() + F::from(3.0).unwrap() / (F::from(2.0).unwrap() * x_sq);
        } else {
            return -(x.cos() / x + x.sin()) / x * x / (-x.cos());
        }
    }

    if x < F::from(5.0).unwrap() {
        let y_n = spherical_yn(n, x);
        return y_n * x / (-x.cos());
    }

    if x > F::from(n * n).unwrap() || x > F::from(1000.0).unwrap() {
        let x_sq = x * x;
        let n_f = F::from(n).unwrap();

        let sign = if n % 2 == 0 { F::one() } else { F::one().neg() };
        let mut factor = sign;

        factor = factor - sign * n_f * (n_f + F::one()) / (F::from(2.0).unwrap() * x_sq);

        if x > F::from(100.0).unwrap() {
            let term2 =
                sign * n_f * (n_f + F::one()) * (n_f * (n_f + F::one()) - F::from(2.0).unwrap())
                    / (F::from(8.0).unwrap() * x_sq * x_sq);
            factor = factor + term2;
        }

        return factor;
    }

    spherical_yn(n, x) * x / (-x.cos())
}
