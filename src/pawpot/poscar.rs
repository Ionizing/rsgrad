// Minimal POSCAR reader (VASP 5+ format).

use anyhow::{bail, Context, Result};
use std::fs;
use std::path::Path;

/// POSCAR structure: cell, element counts, and scaled positions.
#[derive(Debug, Clone)]
pub struct Poscar {
    /// Scale factor
    pub scale: f64,
    /// Lattice vectors [a1, a2, a3], each in Å (after applying scale).
    /// cell[i] is the i-th lattice vector as [x, y, z].
    pub cell: [[f64; 3]; 3],
    /// Element symbols (VASP 5+ format).
    pub symbols: Vec<String>,
    /// Number of atoms per element.
    pub counts: Vec<usize>,
    /// Fractional (scaled) coordinates, shape (natoms, 3).
    pub positions: Vec<[f64; 3]>,
    /// Total number of atoms.
    pub natoms: usize,
}

impl std::str::FromStr for Poscar {
    type Err = anyhow::Error;
    fn from_str(text: &str) -> Result<Self> {
        Self::parse(text)
    }
}

impl Poscar {
    /// Read POSCAR from file.
    pub fn from_file(path: impl AsRef<Path>) -> Result<Self> {
        let text = fs::read_to_string(path.as_ref())
            .with_context(|| format!("reading {}", path.as_ref().display()))?;
        text.parse::<Self>()
    }

    /// Parse POSCAR text.
    fn parse(text: &str) -> Result<Self> {
        let lines: Vec<&str> = text.lines().collect();
        if lines.len() < 8 {
            bail!("POSCAR too short");
        }

        // Line 1: comment (skip)
        // Line 2: scale factor
        let scale: f64 = lines[1]
            .trim()
            .parse()
            .context("parsing scale factor")?;

        // Lines 3-5: lattice vectors
        let mut cell = [[0.0_f64; 3]; 3];
        for (i, row) in cell.iter_mut().enumerate() {
            let vals = parse_floats(lines[2 + i])?;
            if vals.len() < 3 {
                bail!("lattice vector {} has fewer than 3 components", i + 1);
            }
            *row = [vals[0] * scale, vals[1] * scale, vals[2] * scale];
        }

        // Line 6: element symbols (VASP 5+)
        let symbols: Vec<String> = lines[5]
            .split_whitespace()
            .map(|s| s.to_string())
            .collect();

        // Line 7: atom counts
        let counts: Vec<usize> = lines[6]
            .split_whitespace()
            .map(|s| s.parse::<usize>().context("parsing atom count"))
            .collect::<Result<Vec<_>>>()?;

        if symbols.len() != counts.len() {
            bail!("element symbol count != atom count entries");
        }

        let natoms: usize = counts.iter().sum();

        // Line 8: "Direct" or "Cartesian"
        let coord_type = lines[7].trim().to_lowercase();
        let is_direct = coord_type.starts_with('d');

        // Read positions
        let mut positions = Vec::with_capacity(natoms);
        for i in 0..natoms {
            let line_idx = 8 + i;
            if line_idx >= lines.len() {
                bail!("not enough position lines");
            }
            let vals = parse_floats(lines[line_idx])?;
            if vals.len() < 3 {
                bail!("position line {} has fewer than 3 values", i);
            }
            let pos = [vals[0], vals[1], vals[2]];
            if is_direct {
                positions.push(pos);
            } else {
                // Cartesian: convert to fractional using f = B @ cart
                // where B = inv(cell)^T = reciprocal cell (without 2π factor).
                let frac = cart_to_frac(&cell, &pos);
                positions.push(frac);
            }
        }

        Ok(Poscar {
            scale,
            cell,
            symbols,
            counts,
            positions,
            natoms,
        })
    }

    /// Volume of the unit cell (Å³).
    pub fn volume(&self) -> f64 {
        let [a1, a2, a3] = &self.cell;
        // det([[a1],[a2],[a3]])
        a1[0] * (a2[1] * a3[2] - a2[2] * a3[1])
            - a1[1] * (a2[0] * a3[2] - a2[2] * a3[0])
            + a1[2] * (a2[0] * a3[1] - a2[1] * a3[0])
    }

    /// Reciprocal lattice vectors B = (A^{-1})^T, in Å^{-1}.
    /// b_cell[i] is the i-th reciprocal lattice vector.
    pub fn reciprocal_cell(&self) -> [[f64; 3]; 3] {
        let vol = self.volume();
        let [a, b, c] = &self.cell;
        let b0 = cross(b, c);
        let b1 = cross(c, a);
        let b2 = cross(a, b);
        [
            [b0[0] / vol, b0[1] / vol, b0[2] / vol],
            [b1[0] / vol, b1[1] / vol, b1[2] / vol],
            [b2[0] / vol, b2[1] / vol, b2[2] / vol],
        ]
    }

    /// Per-atom element index (0-based).
    pub fn element_idx(&self) -> Vec<usize> {
        self.counts
            .iter()
            .enumerate()
            .flat_map(|(i, &cnt)| std::iter::repeat(i).take(cnt))
            .collect()
    }
}

/// Convert a Cartesian position to fractional coordinates.
/// `f = (A^T)^{-1} @ cart` where A has rows a1, a2, a3.
fn cart_to_frac(cell: &[[f64; 3]; 3], cart: &[f64; 3]) -> [f64; 3] {
    let [a1, a2, a3] = cell;
    let vol = a1[0] * (a2[1] * a3[2] - a2[2] * a3[1])
            - a1[1] * (a2[0] * a3[2] - a2[2] * a3[0])
            + a1[2] * (a2[0] * a3[1] - a2[1] * a3[0]);
    let b0 = cross(a2, a3);
    let b1 = cross(a3, a1);
    let b2 = cross(a1, a2);
    [
        (b0[0] * cart[0] + b0[1] * cart[1] + b0[2] * cart[2]) / vol,
        (b1[0] * cart[0] + b1[1] * cart[1] + b1[2] * cart[2]) / vol,
        (b2[0] * cart[0] + b2[1] * cart[1] + b2[2] * cart[2]) / vol,
    ]
}

fn cross(a: &[f64; 3], b: &[f64; 3]) -> [f64; 3] {
    [
        a[1] * b[2] - a[2] * b[1],
        a[2] * b[0] - a[0] * b[2],
        a[0] * b[1] - a[1] * b[0],
    ]
}

fn parse_floats(line: &str) -> Result<Vec<f64>> {
    // Strip inline comments
    let s = line.split(['!', '#']).next().unwrap_or(line);
    s.split_whitespace()
        .map(|t| t.parse::<f64>().with_context(|| format!("parsing '{t}'")))
        .collect()
}

#[cfg(test)]
mod tests {
    use super::*;

    #[test]
    fn test_poscar_mos2() {
        let p = Poscar::from_file("tests/pawpot/projectors_lreal_false/POSCAR").unwrap();
        assert_eq!(p.symbols, vec!["Mo", "S"]);
        assert_eq!(p.counts, vec![1, 2]);
        assert_eq!(p.natoms, 3);

        let vol = p.volume();
        assert!((vol - 176.43).abs() < 0.5, "volume = {}", vol);

        // Scaled position of Mo at (0, 0, 0)
        assert!(p.positions[0][0].abs() < 1e-6);
        assert!(p.positions[0][1].abs() < 1e-6);
        assert!(p.positions[0][2].abs() < 1e-6);
    }
}
