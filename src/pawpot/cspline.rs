/// Represents the coefficients of the cubic polynomial $S_i(x)$ for the interval $[x_i, x_{i+1}]$.
/// The polynomial is defined as:
/// $S_i(x) = y_i + b_i(x - x_i) + c_i(x - x_i)^2 + d_i(x - x_i)^3$
#[repr(C)]
#[derive(Clone, Copy, Default, Debug)]
struct PolynomialCoefficients {
    /// Coefficient $b_i$
    b: f64,
    /// Coefficient $c_i$ (proportional to $S_i''(x_i)$)
    c: f64,
    /// Coefficient $d_i$
    d: f64,
}

impl AsRef<[f64; 3]> for PolynomialCoefficients {
    fn as_ref(&self) -> &[f64; 3] {
        let ptr = self as *const Self;
        let ptr2 = ptr as *const [f64; 3];
        unsafe { &*ptr2 }
    }
}


impl From<PolynomialCoefficients> for [f64; 3] {
    fn from(coeffs: PolynomialCoefficients) -> Self {
        [coeffs.b, coeffs.c, coeffs.d]
    }
}


/// The main structure for a Cubic Spline interpolator.
/// It stores the original data points and the pre-calculated polynomial coefficients
/// for each segment.
#[derive(Clone)]
pub struct CubicSpline {
    /// The sorted x-coordinates (knots) of the data points.
    /// Required to be strictly monotonic.
    x_knots: Vec<f64>,

    /// The y-coordinates (values) of the data points.
    y_values: Vec<f64>,

    /// The polynomial coefficients for each interval $[x_i, x_{i+1}]$.
    /// The length is $n-1$, where $n$ is the number of knots.
    coefficients: Vec<PolynomialCoefficients>,
}

/// Defines the type of boundary condition for the first or last point.
///
/// * `SecondDerivative(val)`: Specifies the value of the second derivative $S''(x)$ at the boundary.
/// * `FirstDerivative(val)`: Specifies the value of the first derivative $S'(x)$ at the boundary.
#[derive(Clone, Copy)]
pub enum BoundaryCondition {
    SecondDerivative(f64),
    FirstDerivative(f64),
}

impl CubicSpline {
    /// Constructs a **Natural Cubic Spline** from a set of data points $(x, y)$.
    ///
    /// The natural boundary condition sets the second derivative to zero at both endpoints:
    /// $S''(x_0) = 0$ and $S''(x_{n-1}) = 0$.
    ///
    /// # Panics
    /// * If the number of points is less than 3.
    /// * If the `x` coordinates are not monotonic (must be non-decreasing).
    pub fn from_xy_natural(x: &[f64], y: &[f64]) -> Self {
        let bc = [
            BoundaryCondition::SecondDerivative(0.0), // S''(x_0) = 0
            BoundaryCondition::SecondDerivative(0.0), // S''(x_{n-1}) = 0
        ];
        Self::from_xy_with_bc(x, y, &bc)
    }

    /// Constructs a Cubic Spline with the first derivative clamped at the start point
    /// and a natural boundary (second derivative zero) at the end point.
    ///
    /// # Arguments
    /// * `derivative_start`: The desired value of the first derivative $S'(x_0)$ at the start.
    pub fn from_xy_clamped_start(x: &[f64], y: &[f64], derivative_start: f64) -> Self {
        let bc = [
            BoundaryCondition::FirstDerivative(derivative_start),
            BoundaryCondition::SecondDerivative(0.0),
        ];
        Self::from_xy_with_bc(x, y, &bc)
    }

    /// The general constructor for the Cubic Spline, allowing for arbitrary boundary conditions.
    ///
    /// # Arguments
    /// * `x`: The x-coordinates (knots). Must be non-decreasing.
    /// * `y`: The y-coordinates (values).
    /// * `bc`: An array of two boundary conditions: `[start_bc, end_bc]`.
    ///
    /// # Panics
    /// * If the number of points is less than 3.
    /// * If the input x is not monotonic (must be non-decreasing).
    pub fn from_xy_with_bc(x: &[f64], y: &[f64], bc: &[BoundaryCondition; 2]) -> Self {
        let n = x.len(); // Number of knots

        assert_eq!(n, y.len(), "x and y arrays must have the same length");
        assert!(n >= 3, "Cubic Spline requires at least 3 data points");

        // dx[i] = x_{i+1} - x_i
        let delta_x = diff(x);
        // dy[i] = y_{i+1} - y_i
        let delta_y = diff(y);

        assert!(delta_x.iter().all(|&v| v >= 0.0), "Input x coordinates must be monotonic (non-decreasing)");

        // --- Setup for the Tridiagonal System $Ac = r$ (to solve for $c_i = S''(x_i)/2$) ---

        // The system is $Ac = r$, where $c$ is the vector of $c_i$ coefficients.
        // A is a tridiagonal matrix.

        // Lower diagonal $a_i$ (dl), length $n-1$.
        // Note: The loop below implicitly skips $a_0$. We use $a_{i-1}$ in the solve.
        let mut lower_diag = delta_x.clone();
        // The last element of $a$ (corresponding to $i=n-1$) is handled by the boundary condition.
        lower_diag[n-2] = 0.0;

        // Upper diagonal $c_i$ (du), length $n-1$.
        let mut upper_diag = delta_x.clone();
        // The first element of $c$ (corresponding to $i=0$) is handled by the boundary condition.
        upper_diag[0] = 0.0;

        // Main diagonal $b_i$ (d), length $n$.
        let mut main_diag = vec![0.0; n];
        main_diag[0] = 1.0;
        main_diag[n-1] = 1.0;

        // Calculate internal main diagonal elements: $b_i = 2(\Delta x_{i-1} + \Delta x_i)$
        for i in 1 .. n-1 {
            main_diag[i] = 2.0 * (delta_x[i-1] + delta_x[i]);
        }

        // Right-hand side vector $r$, length $n$.
        let mut rhs_vector = vec![0.0; n];
        // Calculate internal RHS elements: $r_i = 3(\frac{\Delta y_i}{\Delta x_i} - \frac{\Delta y_{i-1}}{\Delta x_{i-1}})$
        for i in 1 .. n-1 {
            rhs_vector[i] = 3.0 * (delta_y[i] / delta_x[i] - delta_y[i-1] / delta_x[i-1]);
        }

        // --- Apply Boundary Conditions (BC) ---

        // BC at $x_0$ (i=0)
        match bc[0] {
            // Natural Spline or Second Derivative Clamped: $S''(x_0) = 2c_0 = val$.
            BoundaryCondition::SecondDerivative(val) => {
                // Simplified system for $2c_0 = val$:
                // $1 \cdot c_0 + 0 \cdot c_1 = val/2.0$
                // However, the Thomas algorithm expects $b_0 c_0 + c_0 c_1 = r_0$.
                // For $c_0 = val/2.0$, we simplify the first row to:
                // $2 c_0 + 0 c_1 = val$. This means $b_0=2$, $r_0=val$.
                main_diag[0] = 2.0;
                rhs_vector[0] = val;
                // $c_0$ is zero for $i=0$.
                upper_diag[0] = 0.0;
            },
            // First Derivative Clamped: $S'(x_0) = b_0 = val$.
            // $2 \Delta x_0 c_0 + \Delta x_0 c_1 = 3 (\frac{\Delta y_0}{\Delta x_0} - val)$
            BoundaryCondition::FirstDerivative(val) => {
                main_diag[0] = 2.0 * delta_x[0]; // $b_0$
                upper_diag[0] = delta_x[0];      // $c_0$
                rhs_vector[0] = 3.0 * (delta_y[0] / delta_x[0] - val); // $r_0$
            }
        }

        // BC at $x_{n-1}$ (i=n-1)
        match bc[1] {
            // Natural Spline or Second Derivative Clamped: $S''(x_{n-1}) = 2c_{n-1} = val$.
            BoundaryCondition::SecondDerivative(val) => {
                // Simplified system for $2c_{n-1} = val$:
                // $0 c_{n-2} + 2 c_{n-1} = val$
                main_diag[n-1] = 2.0;
                rhs_vector[n-1] = val;
                // $a_{n-2}$ is zero.
                lower_diag[n-2] = 0.0;
            },
            // First Derivative Clamped: $S'(x_{n-1}) = b_{n-1} = val$.
            // $\Delta x_{n-2} c_{n-2} + 2 \Delta x_{n-2} c_{n-1} = 3 (val - \frac{\Delta y_{n-2}}{\Delta x_{n-2}})$
            BoundaryCondition::FirstDerivative(val) => {
                main_diag[n-1] = 2.0 * delta_x[n-2]; // $b_{n-1}$
                lower_diag[n-2] = delta_x[n-2];      // $a_{n-2}$
                rhs_vector[n-1] = 3.0 * (val - delta_y[n-2] / delta_x[n-2]); // $r_{n-1}$
            }
        }

        // --- Solve the System ---

        // Solve for the $c$ coefficients using the Thomas Algorithm.
        let c_coefficients = solve_tridiagonal(&lower_diag, &main_diag, &upper_diag, &rhs_vector);

        // --- Calculate $b$ and $d$ coefficients ---

        let mut polynomials = vec![PolynomialCoefficients::default(); n-1];
        for i in 0 .. n-1 {
            let dx_i = delta_x[i];
            let dy_i = delta_y[i];
            let c_i = c_coefficients[i];
            let c_i_plus_1 = c_coefficients[i+1];

            // $c_i$ is already calculated (proportional to $S''(x_i)$).
            polynomials[i].c = c_i;

            // $b_i = \frac{\Delta y_i}{\Delta x_i} - \frac{\Delta x_i}{3}(2c_i + c_{i+1})$
            polynomials[i].b = dy_i / dx_i - dx_i / 3.0 * (2.0 * c_i + c_i_plus_1);

            // $d_i = \frac{c_{i+1} - c_i}{3 \Delta x_i}$
            polynomials[i].d = (c_i_plus_1 - c_i) / (3.0 * dx_i);
        }

        Self {
            x_knots: x.to_owned(),
            y_values: y.to_owned(),
            coefficients: polynomials,
        }
    }

    /// Evaluates the cubic spline function $S(x_0)$ at a given point $x_0$.
    ///
    /// # Arguments
    /// * `x0`: The point at which to evaluate the spline.
    ///
    /// # Returns
    /// * `Some(f64)`: The interpolated value $S(x_0)$ if $x_0$ is within $[x_0, x_{n-1}]$.
    /// * `None`: If $x_0$ is outside the range $[x_0, x_{n-1}]$.
    pub fn eval(&self, x0: f64) -> Option<f64> {
        // Handle the exact start point
        if x0 == self.x_knots[0] {
            return Some(self.y_values[0]);
        }

        // Find the index $i$ such that $x_i \le x_0 < x_{i+1}$.
        // partition_point returns the first index $i0$ for which the predicate is false,
        // i.e., the first element not less than $x0$.
        let i0 = self.x_knots.as_slice().partition_point(|v| *v < x0);

        // x0 is out of range:
        // Case 1: $x0 < x_0$. i0 will be 0.
        // Case 2: $x0 > x_{n-1}$. i0 will be $n$.
        if 0 == i0 || i0 >= self.x_knots.len() {
            return None;
        }

        // The segment index is $i = i0 - 1$.
        let segment_index = i0 - 1;

        if x0 == self.x_knots[segment_index + 1] {
            return Some(self.y_values[segment_index + 1]);
        }

        // Local coordinate: $\xi = x_0 - x_i$
        let xi = x0 - self.x_knots[segment_index];

        // $y_i$
        let yi = self.y_values[segment_index];

        // Polynomial coefficients for segment $i$
        let pi = self.coefficients[segment_index];

        // Evaluate $S_i(x_0) = y_i + b_i\xi + c_i\xi^2 + d_i\xi^3$
        // This uses Horner's method for efficiency and stability:
        // $y_i + \xi (b_i + \xi (c_i + d_i \xi))$
        Some( yi + xi * (pi.b + xi * (pi.c + pi.d * xi) ) )
    }

    pub fn eval_vec(&self, xi: &[f64]) -> Option<Vec<f64>> {
        xi.iter().map(|&x| self.eval(x)).collect()
    }

    pub fn dump_coefficients(&self) -> Vec<[f64; 3]> {
        self.coefficients.iter()
            .map(|&x| x.into())
            .collect()
    }
}

/// Calculates the difference between consecutive elements: $z_i = x_{i+1} - x_i$.
fn diff(x: &[f64]) -> Vec<f64> {
    std::iter::zip(x.iter().cloned(), x[1..].iter().cloned())
        .map(|(x1, x2)| x2 - x1)
        .collect()
}


/// Solves a tridiagonal linear system $Ax = r$ using the Thomas Algorithm (TDMA).
fn solve_tridiagonal(lower_diag: &[f64], main_diag: &[f64], upper_diag: &[f64], rhs_vector: &[f64]) -> Vec<f64> {
    let n = main_diag.len(); // Matrix size $n \times n$

    // 1. Input Validation
    assert!(n >= 2, "Matrix size must be at least 2x2");
    assert_eq!(lower_diag.len() + 1, n, "Lower diagonal length mismatch (expected n-1)");
    assert_eq!(upper_diag.len() + 1, n, "Upper diagonal length mismatch (expected n-1)");
    assert_eq!(rhs_vector.len(), n, "RHS vector length mismatch (expected n)");

    // `x` will initially store the modified RHS ($r'$ or $d'$) and later the solution vector.
    let mut x = rhs_vector.to_owned();

    // `c_prime` stores the modified upper diagonal coefficients ($c'$).
    let mut c_prime = vec![0.0; n];

    // --- Forward Elimination Phase (Tridiagonal LU Decomposition) ---

    // Row $i=0$ (Handle separately): $c'_0 = c_0 / b_0$, $r'_0 = r_0 / b_0$
    if main_diag[0] == 0.0 {
        panic!("Singular matrix: zero pivot at row 0");
    }
    // Only calculate $c'_0$ if $i < n-1$.
    if n > 1 {
        c_prime[0] = upper_diag[0] / main_diag[0];
    }
    x[0] /= main_diag[0];

    // Rows $i=1$ to $n-1$.
    for i in 1..n {
        // Calculate the denominator (pivot element): $b'_i = b_i - a_{i-1} c'_{i-1}$
        let pivot_denom = main_diag[i] - lower_diag[i - 1] * c_prime[i - 1];

        if pivot_denom == 0.0 {
            panic!("Singular matrix: zero pivot at row {}", i);
        }

        // Calculate new upper diagonal coefficient $c'_i$.
        // This is only needed for $i < n-1$.
        if i < n - 1 {
            // $c'_i = c_i / b'_i$
            c_prime[i] = upper_diag[i] / pivot_denom;
        }

        // Calculate new RHS value $r'_i$.
        // $r'_i = (r_i - a_{i-1} r'_{i-1}) / b'_i$
        x[i] = (x[i] - lower_diag[i - 1] * x[i - 1]) / pivot_denom;
    }

    // --- Backward Substitution Phase ---

    // The last element $x_{n-1}$ is already solved: $x_{n-1} = r'_{n-1}$.
    // Solve for $x_i$ by iterating backwards from $n-2$ to 0.
    // Formula: $x_i = r'_i - c'_i x_{i+1}$
    for i in (0..n - 1).rev() {
        x[i] -= c_prime[i] * x[i + 1];
    }

    x
}


#[cfg(test)]
mod tests {
    use std::collections::HashMap;
    use std::sync::OnceLock;
    use approx::assert_abs_diff_eq;
    use serde::Deserialize;
    use super::*;


    #[derive(Debug, Deserialize)]
    struct SplineTestRef {
        x_knots: Vec<f64>,
        y_values: Vec<f64>,
        coefficients: Vec<[f64; 3]>,
        x_interpolated: Vec<f64>,
        y_interpolated: Vec<f64>,
    }


    /// Generates a vector of `num` evenly spaced samples over the interval [`start`, `stop`].
    fn linspace(start: f64, stop: f64, num: usize) -> Vec<f64> {
        match num {
            0 => return Vec::new(),
            1 => return vec![start],
            _ => {}
        }

        let step = (stop - start) / ((num - 1) as f64);
        let mut result = Vec::with_capacity(num);

        for i in 0..num {
            let value = start + (i as f64) * step;
            result.push(value);
        }

        // Ensure the last element is exactly 'stop' to prevent floating-point drift.
        if num > 0 {
            *result.last_mut().unwrap() = stop;
        }

        result
    }

    fn get_testset_from_json() -> &'static HashMap<String, SplineTestRef>
    {
        static CUBIC_SPLINE_TESTS: OnceLock<HashMap<String, SplineTestRef>> = OnceLock::new();
        const TEST_FNAME: &str = concat!(env!("CARGO_MANIFEST_DIR"), "/tests/pawpot/cspline_ref.json");
        &CUBIC_SPLINE_TESTS.get_or_init(|| {
            let content = std::fs::read_to_string(TEST_FNAME)
                .expect(&format!("Could not read reference file: {:?}", TEST_FNAME));
            let data: HashMap<String, SplineTestRef> = serde_json::from_str(&content)
                .expect("Failed to deserialize JSON data into HashMap");
            data
        })
    }


    #[test]
    fn test_solve_tridiagonal_simple() {
        // System:
        // 2x + 1y = 8
        // 1x + 3y + 1z = 12.5
        //      1y + 2z = 9
        // Solution should be x=3, y=2, z=3.5
        let lower_diag = vec![1.0, 1.0];
        let main_diag = vec![2.0, 3.0, 2.0];
        let upper_diag = vec![1.0, 1.0];
        let rhs_vector = vec![8.0, 12.5, 9.0];
        let solution = solve_tridiagonal(&lower_diag, &main_diag, &upper_diag, &rhs_vector);
        let expected = vec![3.0f64, 2.0, 3.5];

        assert_abs_diff_eq!(solution.as_slice(), expected.as_slice(), epsilon=f64::EPSILON * 2.0);
    }

    #[test]
    fn test_cspline_natural() {
        let x = vec![1.0, 2.0, 3.0, 4.0];
        let y = vec![0.5, 1.2, 3.4, 2.5];

        let cs = CubicSpline::from_xy_natural(&x, &y);

        // Exact knots
        assert_eq!(cs.eval(1.0).unwrap(), 0.5);
        assert_eq!(cs.eval(2.0).unwrap(), 1.2);
        assert_eq!(cs.eval(3.0).unwrap(), 3.4);
        assert_eq!(cs.eval(4.0).unwrap(), 2.5);

        // Compare coefficients directly
        let coefficients: Vec<[f64; 3]> = vec![
            [ 9.33333333e-02,  1.11022302e-16,  6.06666667e-01],
            [ 1.91333333e+00,  1.82000000e+00, -1.53333333e+00],
            [ 9.53333333e-01, -2.78000000e+00,  9.26666667e-01],
        ];
        let cs_flat: Vec<f64> = cs.dump_coefficients().into_iter().flatten().collect();
        let coeff_flat: Vec<f64> = coefficients.into_iter().flatten().collect();
        assert_abs_diff_eq!(cs_flat.as_slice(), coeff_flat.as_slice(), epsilon=1E-8);

    }

    #[test]
    fn test_cspline_natural_with_testset() {
        fn helper(casename: &str) {
            let simple = &get_testset_from_json()[casename];
            let x_knots = &simple.x_knots;
            let y_values = &simple.y_values;
            let cs = CubicSpline::from_xy_natural(x_knots, y_values);

            // Test exact knots.
            let yi = cs.eval_vec(x_knots).unwrap();
            assert_eq!(y_values, &yi);

            // Check Coefficients
            let coefficients = cs.dump_coefficients();
            let cs_flat: Vec<f64> = coefficients.into_iter().flatten().collect();
            let coeff_flat: Vec<f64> = simple.coefficients.iter().flatten().copied().collect();
            assert_abs_diff_eq!(cs_flat.as_slice(), coeff_flat.as_slice(), epsilon=1E-12);

            let xi = &simple.x_interpolated;
            let yi = cs.eval_vec(xi).unwrap();
            assert_abs_diff_eq!(yi.as_slice(), simple.y_interpolated.as_slice(), epsilon=1E-12);
        }

        helper("cubic_spline_natural_simple");
        helper("cubic_spline_natural_many");
    }
}
