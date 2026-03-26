// All-electron wavefunction reconstruction (aewfc).

use std::f64::consts::{PI, SQRT_2, TAU};
use std::path::Path;

use super::npy::{write_npy_c128, NpzWriter};

use anyhow::Result;
use num::complex::Complex;

use super::gvec::{gvectors, gvectors_gamma};
use super::nonlq::Nonlq;
use super::pawpotcar::PawPotcar;
use super::poscar::Poscar;
use super::sph_harm::sph_r;
use super::spherical_bessel::spherical_jn;
use super::wavecar::Wavecar;

type C128 = Complex<f64>;
/// Per-type SBT projector array: [ntypes][lmmax][n_ae_gvec].
type SbtProjArray = Vec<Vec<Vec<f64>>>;

const TPI: f64 = TAU;
const AUTOA: f64 = 0.529177249;
const RYTOEV: f64 = 13.605826;

/// All-electron wavefunction reconstruction for a single k-point.
pub struct VaspAeWfc {
    wavecar: Wavecar,
    ikpt: usize,
    pub lgam: bool,

    pub nonlq: Nonlq,
    ps_gvec: Vec<[i32; 3]>,

    pub aegrid: [usize; 3],
    ae_gvec: Vec<[i32; 3]>,
    #[allow(dead_code)]
    ae_glen: Vec<f64>,

    element_idx: Vec<usize>,
    lmmax_per_type: Vec<usize>,

    qij_per_type: Vec<Vec<Vec<f64>>>,

    crexp_ae: Vec<Vec<C128>>,
    cqfak_ae: Vec<Vec<C128>>,

    q_ae_core: Vec<Vec<Vec<f64>>>,
    q_ps_core: Vec<Vec<Vec<f64>>>,
}

impl VaspAeWfc {
    /// Build the AE-WFC reconstruction object.
    pub fn new(
        wavecar: Wavecar,
        poscar: &Poscar,
        pawpot: &PawPotcar,
        ikpt: usize,
        aecut: f64,
    ) -> Result<Self> {
        let pscut = wavecar.header.encut;
        let kvec = wavecar.kpoints[ikpt - 1].kvec;
        let bcell = wavecar.bcell();
        let volume = poscar.volume().abs();
        let aecut = if aecut < 0.0 { -aecut * pscut } else { aecut.max(pscut) };

        let ngrid_ps = compute_ngrid_from_cell(&poscar.cell, pscut);
        let nplw_stored = wavecar.kpoints[ikpt - 1].nplw;
        let nplw_full = gvectors(&bcell, pscut, &kvec, &ngrid_ps).0.len();
        let nplw_gamma = gvectors_gamma(&bcell, pscut, &kvec, &ngrid_ps).0.len();
        let lgam = nplw_stored == nplw_gamma && nplw_stored != nplw_full;

        let nonlq = Nonlq::new(poscar, pawpot, pscut, &kvec, &bcell, lgam);
        let ps_gvec = nonlq.gvec.clone();

        let aegrid = compute_aegrid(aecut, pscut, &wavecar.header.ngrid);
        let (ae_gvec, ae_gk, ae_glen) = if lgam {
            gvectors_gamma(&bcell, aecut, &kvec, &aegrid)
        } else {
            gvectors(&bcell, aecut, &kvec, &aegrid)
        };

        let natoms = poscar.natoms;
        let element_idx = poscar.element_idx();
        let ntypes = pawpot.len();

        let lmmax_per_type: Vec<usize> = (0..ntypes)
            .map(|it| pawpot.get(it).unwrap().lmmax())
            .collect();
        let qij_per_type: Vec<Vec<Vec<f64>>> = (0..ntypes)
            .map(|it| pawpot.get(it).unwrap().get_qij())
            .collect();

        let n_ae = ae_gvec.len();
        let scaled_pos = &poscar.positions;
        let mut crexp_ae = vec![vec![C128::new(0.0, 0.0); natoms]; n_ae];
        for (ig, g) in ae_gvec.iter().enumerate() {
            for (iatom, pos) in scaled_pos.iter().enumerate() {
                let dot = g[0] as f64 * pos[0] + g[1] as f64 * pos[1] + g[2] as f64 * pos[2];
                crexp_ae[ig][iatom] = Complex::new(0.0, -TPI * dot).exp();
            }
        }

        let cqfak_ae: Vec<Vec<C128>> = (0..ntypes).map(|it| {
            let pp = pawpot.get(it).unwrap();
            pp.ls().iter()
                .flat_map(|&l| {
                    let il_inv = C128::new(0.0, 1.0).powi(-(l as i32));
                    vec![il_inv; 2 * l + 1]
                })
                .collect()
        }).collect();

        let lmax_global = (0..ntypes)
            .map(|it| pawpot.get(it).unwrap().ls().iter().copied().max().unwrap_or(0))
            .max()
            .unwrap_or(0);
        let ylm_ae: Vec<Vec<Vec<f64>>> = (0..=lmax_global)
            .map(|l| sph_r(&ae_gk, l))
            .collect();

        let (q_ae_core, q_ps_core) = build_sbt_projectors(
            pawpot, ntypes, &ae_glen, &ylm_ae, volume,
        )?;

        Ok(VaspAeWfc {
            wavecar,
            ikpt,
            lgam,
            nonlq,
            ps_gvec,
            aegrid,
            ae_gvec,
            ae_glen,
            element_idx,
            lmmax_per_type,
            qij_per_type,
            crexp_ae,
            cqfak_ae,
            q_ae_core,
            q_ps_core,
        })
    }

    /// Compute the AE norm: `Σ|C_G|² + β† Q β`.
    pub fn get_ae_norm(&mut self, ispin: usize, iband: usize) -> Result<f64> {
        let coeffs = self.wavecar.read_band_coeff(ispin, self.ikpt, iband)?;
        let beta = self.nonlq.proj(&coeffs);

        let ps_norm: f64 = coeffs.iter().map(|c| c.norm_sqr()).sum();

        let natoms = self.element_idx.len();
        let mut qij_sum = C128::new(0.0, 0.0);
        let mut nproj = 0;
        for iatom in 0..natoms {
            let ntype = self.element_idx[iatom];
            let lmmax = self.lmmax_per_type[ntype];
            let qij = &self.qij_per_type[ntype];
            for i in 0..lmmax {
                for j in 0..lmmax {
                    qij_sum += beta[nproj + i].conj() * qij[i][j] * beta[nproj + j];
                }
            }
            nproj += lmmax;
        }

        Ok(ps_norm + qij_sum.re)
    }

    /// Reconstruct the AE wavefunction in real space.
    pub fn get_ae_wfc(&mut self, ispin: usize, iband: usize, norm: bool) -> Result<Vec<C128>> {
        let coeffs_raw = self.wavecar.read_band_coeff(ispin, self.ikpt, iband)?;
        let beta = self.nonlq.proj(&coeffs_raw);

        let mut coeffs = coeffs_raw;
        if self.lgam {
            for c in coeffs[1..].iter_mut() {
                *c /= SQRT_2;
            }
        }

        let n_ae = self.ae_gvec.len();
        let natoms = self.element_idx.len();
        let [nx, ny, nz] = self.aegrid;

        let mut sg_ae = vec![C128::new(0.0, 0.0); n_ae];
        let mut sg_ps = vec![C128::new(0.0, 0.0); n_ae];

        let mut nproj = 0;
        for iatom in 0..natoms {
            let ntype = self.element_idx[iatom];
            let lmmax = self.lmmax_per_type[ntype];
            let cqfak = &self.cqfak_ae[ntype];
            let qae = &self.q_ae_core[ntype];
            let qps = &self.q_ps_core[ntype];

            for ig in 0..n_ae {
                let phase = self.crexp_ae[ig][iatom];
                let mut sum_ae = C128::new(0.0, 0.0);
                let mut sum_ps = C128::new(0.0, 0.0);
                for ilm in 0..lmmax {
                    let b = beta[nproj + ilm];
                    let il = cqfak[ilm];
                    sum_ae += il * b * qae[ilm][ig];
                    sum_ps += il * b * qps[ilm][ig];
                }
                sg_ae[ig] += phase * sum_ae;
                sg_ps[ig] += phase * sum_ps;
            }

            nproj += lmmax;
        }

        let mut phi_ae = vec![C128::new(0.0, 0.0); nx * ny * nz];

        let flat_idx = |gx: i32, gy: i32, gz: i32| -> usize {
            let ix = gx.rem_euclid(nx as i32) as usize;
            let iy = gy.rem_euclid(ny as i32) as usize;
            let iz = gz.rem_euclid(nz as i32) as usize;
            ix * ny * nz + iy * nz + iz
        };

        for (ig, g) in self.ps_gvec.iter().enumerate() {
            let idx = flat_idx(g[0], g[1], g[2]);
            phi_ae[idx] += coeffs[ig];
        }

        for (ig, g) in self.ae_gvec.iter().enumerate() {
            let idx = flat_idx(g[0], g[1], g[2]);
            phi_ae[idx] += sg_ae[ig] - sg_ps[ig];
        }

        if self.lgam {
            for (ig, g) in self.ps_gvec.iter().enumerate().skip(1) {
                let idx = flat_idx(-g[0], -g[1], -g[2]);
                phi_ae[idx] += coeffs[ig].conj();
            }
            for (ig, g) in self.ae_gvec.iter().enumerate().skip(1) {
                let idx = flat_idx(-g[0], -g[1], -g[2]);
                phi_ae[idx] += (sg_ae[ig] - sg_ps[ig]).conj();
            }
        }

        ifft3d(&mut phi_ae, nx, ny, nz);

        let total = (nx * ny * nz) as f64;
        if norm {
            let fac = total.sqrt();
            for v in phi_ae.iter_mut() {
                *v /= fac;
            }
        } else {
            for v in phi_ae.iter_mut() {
                *v /= total;
            }
        }

        Ok(phi_ae)
    }

    /// Debug helper: return the raw (sg_ae, sg_ps) G-space vectors for a given band.
    /// These are the on-site AE and PS correction amplitudes over the AE G-sphere,
    /// before the IFFT into real space.
    pub fn compute_sg_ae_ps(
        &mut self,
        ispin: usize,
        iband: usize,
    ) -> Result<(Vec<C128>, Vec<C128>)> {
        let coeffs_raw = self.wavecar.read_band_coeff(ispin, self.ikpt, iband)?;
        let beta = self.nonlq.proj(&coeffs_raw);

        let n_ae = self.ae_gvec.len();
        let natoms = self.element_idx.len();

        let mut sg_ae = vec![C128::new(0.0, 0.0); n_ae];
        let mut sg_ps = vec![C128::new(0.0, 0.0); n_ae];

        let mut nproj = 0;
        for iatom in 0..natoms {
            let ntype = self.element_idx[iatom];
            let lmmax = self.lmmax_per_type[ntype];
            let cqfak = &self.cqfak_ae[ntype];
            let qae = &self.q_ae_core[ntype];
            let qps = &self.q_ps_core[ntype];

            for ig in 0..n_ae {
                let phase = self.crexp_ae[ig][iatom];
                let mut sum_ae = C128::new(0.0, 0.0);
                let mut sum_ps = C128::new(0.0, 0.0);
                for ilm in 0..lmmax {
                    let b = beta[nproj + ilm];
                    let il = cqfak[ilm];
                    sum_ae += il * b * qae[ilm][ig];
                    sum_ps += il * b * qps[ilm][ig];
                }
                sg_ae[ig] += phase * sum_ae;
                sg_ps[ig] += phase * sum_ps;
            }

            nproj += lmmax;
        }

        Ok((sg_ae, sg_ps))
    }

    /// Write the AE wavefunction density |ψ|² to a VASP CHGCAR-format file.
    pub fn write_ae_density(
        &mut self,
        poscar: &super::poscar::Poscar,
        ispin: usize,
        iband: usize,
        path: impl AsRef<std::path::Path>,
    ) -> Result<()> {
        use std::io::Write;

        let phi = self.get_ae_wfc(ispin, iband, true)?;
        let [nx, ny, nz] = self.aegrid;
        let ae_factor = (nx * ny * nz) as f64;

        let mut f = std::io::BufWriter::new(std::fs::File::create(path.as_ref())?);

        writeln!(f, "AE density  ispin={ispin} iband={iband}")?;
        writeln!(f, "1.0")?;
        for v in &poscar.cell {
            writeln!(f, "  {:18.10}  {:18.10}  {:18.10}", v[0], v[1], v[2])?;
        }
        write!(f, " ")?;
        for s in &poscar.symbols {
            write!(f, " {s}")?;
        }
        writeln!(f)?;
        write!(f, " ")?;
        for c in &poscar.counts {
            write!(f, " {c}")?;
        }
        writeln!(f)?;
        writeln!(f, "Direct")?;
        for p in &poscar.positions {
            writeln!(f, "  {:18.10}  {:18.10}  {:18.10}", p[0], p[1], p[2])?;
        }
        writeln!(f)?;

        writeln!(f, "  {nx}  {ny}  {nz}")?;

        let mut count = 0usize;
        for ix in 0..nx {
            for iy in 0..ny {
                for iz in 0..nz {
                    let idx = ix * ny * nz + iy * nz + iz;
                    let rho = phi[idx].norm_sqr() * ae_factor;
                    write!(f, "  {rho:18.10E}")?;
                    count += 1;
                    if count % 5 == 0 {
                        writeln!(f)?;
                    }
                }
            }
        }
        if count % 5 != 0 {
            writeln!(f)?;
        }

        Ok(())
    }

    /// Write the AE wavefunction to a `.npy` file as a complex128 array of
    /// shape `[nx, ny, nz]`.
    pub fn write_ae_wfc_npy(
        &mut self,
        ispin: usize,
        iband: usize,
        path: impl AsRef<Path>,
    ) -> Result<()> {
        let wfc = self.get_ae_wfc(ispin, iband, false)?;
        let [nx, ny, nz] = self.aegrid;
        let shape = [nx as u64, ny as u64, nz as u64];

        write_npy_c128(path.as_ref(), &wfc, &shape)?;
        Ok(())
    }

    /// Write the AE wavefunction to a `.npz` file.
    pub fn write_ae_wfc_npz(
        &mut self,
        ispin: usize,
        iband: usize,
        path: impl AsRef<Path>,
    ) -> Result<()> {
        let wfc = self.get_ae_wfc(ispin, iband, true)?;
        let [nx, ny, nz] = self.aegrid;
        let ae_factor = (nx * ny * nz) as f64;
        let density: Vec<f64> = wfc.iter().map(|c| c.norm_sqr() * ae_factor).collect();
        let shape = [nx as u64, ny as u64, nz as u64];

        let mut npz = NpzWriter::create(path.as_ref())?;
        npz.write_c128("wfc", &wfc, &shape)?;
        npz.write_f64("density", &density, &shape)?;
        npz.finish()?;
        Ok(())
    }
}

// ---------------------------------------------------------------------------
// SBT projector builder
// ---------------------------------------------------------------------------

fn build_sbt_projectors(
    pawpot: &PawPotcar,
    ntypes: usize,
    ae_glen: &[f64],
    ylm_ae: &[Vec<Vec<f64>>],
    volume: f64,
) -> Result<(SbtProjArray, SbtProjArray)> {
    let n_ae = ae_glen.len();
    let inv_sqrt_vol = 1.0 / volume.sqrt();

    let unique_glen = unique_sorted(ae_glen, 1e-10);
    let glen_to_idx: Vec<usize> = ae_glen
        .iter()
        .map(|&g| upper_bound(&unique_glen, g, 1e-10))
        .collect();

    let mut q_ae_out = Vec::with_capacity(ntypes);
    let mut q_ps_out = Vec::with_capacity(ntypes);

    for itype in 0..ntypes {
        let pp = pawpot.get(itype).unwrap();
        let radgrid = pp.radgrid();
        let simp_w = pp.simp_weights();
        let rmax = pp.rmax();
        let lmmax = pp.lmmax();

        let mut t_ae = vec![vec![0.0f64; n_ae]; lmmax];
        let mut t_ps = vec![vec![0.0f64; n_ae]; lmmax];

        let mut ilm = 0usize;
        for (ii, &l) in pp.ls().iter().enumerate() {
            let tlp1 = 2 * l + 1;
            let ae_wave = &pp.aewaves()[ii];
            let ps_wave = &pp.pswaves()[ii];

            let sbt_ae_uniq: Vec<f64> = unique_glen
                .iter()
                .map(|&q| sbt_direct(q, l, radgrid, simp_w, ae_wave, rmax))
                .collect();
            let sbt_ps_uniq: Vec<f64> = unique_glen
                .iter()
                .map(|&q| sbt_direct(q, l, radgrid, simp_w, ps_wave, rmax))
                .collect();

            for (ig, &iu) in glen_to_idx.iter().enumerate().take(n_ae) {
                let gae_val = sbt_ae_uniq[iu] * inv_sqrt_vol;
                let gps_val = sbt_ps_uniq[iu] * inv_sqrt_vol;
                let ylm_vals = &ylm_ae[l][ig];
                for (im, &y) in ylm_vals.iter().enumerate() {
                    t_ae[ilm + im][ig] = gae_val * y;
                    t_ps[ilm + im][ig] = gps_val * y;
                }
            }

            ilm += tlp1;
        }

        q_ae_out.push(t_ae);
        q_ps_out.push(t_ps);
    }

    Ok((q_ae_out, q_ps_out))
}

// ---------------------------------------------------------------------------
// Spherical Bessel Transform (direct integration)
// ---------------------------------------------------------------------------

fn sbt_direct(
    q: f64,
    l: usize,
    radgrid: &[f64],
    simp_weights: &[f64],
    phi: &[f64],
    rmax: f64,
) -> f64 {
    let mut sum = 0.0;
    for i in 0..radgrid.len() {
        if radgrid[i] >= rmax {
            break;
        }
        sum += simp_weights[i] * radgrid[i] * phi[i] * spherical_jn(l as i32, q * radgrid[i]);
    }
    4.0 * PI * sum
}

// ---------------------------------------------------------------------------
// 3D IFFT
// ---------------------------------------------------------------------------

fn ifft3d(data: &mut [C128], nx: usize, ny: usize, nz: usize) {
    use rustfft::FftPlanner;

    let mut planner = FftPlanner::<f64>::new();

    {
        let ifft_z = planner.plan_fft_inverse(nz);
        let mut scratch = vec![C128::new(0.0, 0.0); ifft_z.get_inplace_scratch_len()];
        for ix in 0..nx {
            for iy in 0..ny {
                let start = ix * ny * nz + iy * nz;
                ifft_z.process_with_scratch(&mut data[start..start + nz], &mut scratch);
            }
        }
    }

    {
        let ifft_y = planner.plan_fft_inverse(ny);
        let mut scratch = vec![C128::new(0.0, 0.0); ifft_y.get_inplace_scratch_len()];
        let mut tmp = vec![C128::new(0.0, 0.0); ny];
        for ix in 0..nx {
            for iz in 0..nz {
                for iy in 0..ny {
                    tmp[iy] = data[ix * ny * nz + iy * nz + iz];
                }
                ifft_y.process_with_scratch(&mut tmp, &mut scratch);
                for iy in 0..ny {
                    data[ix * ny * nz + iy * nz + iz] = tmp[iy];
                }
            }
        }
    }

    {
        let ifft_x = planner.plan_fft_inverse(nx);
        let mut scratch = vec![C128::new(0.0, 0.0); ifft_x.get_inplace_scratch_len()];
        let mut tmp = vec![C128::new(0.0, 0.0); nx];
        for iy in 0..ny {
            for iz in 0..nz {
                for ix in 0..nx {
                    tmp[ix] = data[ix * ny * nz + iy * nz + iz];
                }
                ifft_x.process_with_scratch(&mut tmp, &mut scratch);
                for ix in 0..nx {
                    data[ix * ny * nz + iy * nz + iz] = tmp[ix];
                }
            }
        }
    }
}

// ---------------------------------------------------------------------------
// Helpers
// ---------------------------------------------------------------------------

fn compute_aegrid(aecut: f64, pscut: f64, ngrid: &[usize; 3]) -> [usize; 3] {
    let fac = (aecut / pscut).sqrt();
    [
        (fac * ngrid[0] as f64 + 1.0) as usize,
        (fac * ngrid[1] as f64 + 1.0) as usize,
        (fac * ngrid[2] as f64 + 1.0) as usize,
    ]
}

fn compute_ngrid_from_cell(cell: &[[f64; 3]; 3], encut: f64) -> [usize; 3] {
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

fn unique_sorted(values: &[f64], tol: f64) -> Vec<f64> {
    let mut sorted = values.to_vec();
    sorted.sort_by(|a, b| a.partial_cmp(b).unwrap());
    let mut unique = Vec::with_capacity(sorted.len());
    let mut prev = f64::NEG_INFINITY;
    for &v in &sorted {
        if (v - prev).abs() > tol {
            unique.push(v);
            prev = v;
        }
    }
    unique
}

fn upper_bound(unique: &[f64], val: f64, tol: f64) -> usize {
    match unique.binary_search_by(|&u| {
        if (u - val).abs() <= tol {
            std::cmp::Ordering::Equal
        } else {
            u.partial_cmp(&val).unwrap()
        }
    }) {
        Ok(i) => i,
        Err(i) => {
            let lo = if i > 0 { i - 1 } else { 0 };
            let hi = if i < unique.len() { i } else { unique.len() - 1 };
            if (unique[lo] - val).abs() <= (unique[hi] - val).abs() { lo } else { hi }
        }
    }
}
