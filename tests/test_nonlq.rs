// Integration test for nonlq reciprocal-space projectors.
// Compares Rust Nonlq against reference values from beta.npy (computed by paw.py).

use std::convert::TryInto;

use rsgrad::pawpot::{Nonlq, PawPotcar, PawPoscar, PawWavecar};
use rsgrad::pawpot::VaspAeWfc;
use num::complex::Complex;

type C128 = Complex<f64>;

/// Minimal .npy reader for complex128 arrays.
/// Only handles 1-D and 2-D arrays with dtype '<c16' (little-endian complex128).
fn read_npy_complex128(path: &str) -> Vec<C128> {
    let bytes = std::fs::read(path).expect("read npy file");
    // Magic: \x93NUMPY
    assert_eq!(&bytes[0..6], b"\x93NUMPY", "bad npy magic");
    let major = bytes[6];
    let header_len = if major == 1 {
        u16::from_le_bytes([bytes[8], bytes[9]]) as usize
    } else {
        u32::from_le_bytes([bytes[8], bytes[9], bytes[10], bytes[11]]) as usize
    };
    let header_offset = if major == 1 { 10 } else { 12 };
    let _header = std::str::from_utf8(&bytes[header_offset..header_offset + header_len])
        .expect("utf8 header");
    let data_start = header_offset + header_len;
    let raw = &bytes[data_start..];
    let n = raw.len() / 16; // 16 bytes per complex128
    (0..n)
        .map(|i| {
            let re = f64::from_le_bytes(raw[16 * i..16 * i + 8].try_into().unwrap());
            let im = f64::from_le_bytes(raw[16 * i + 8..16 * i + 16].try_into().unwrap());
            C128::new(re, im)
        })
        .collect()
}

#[test]
fn test_nonlq_proj_vs_reference() {
    let dir = "tests/pawpot/projectors_lreal_false";

    // Load structures
    let poscar = PawPoscar::from_file(format!("{dir}/POSCAR")).unwrap();
    let pawpot = PawPotcar::from_file(&format!("{dir}/POTCAR")).unwrap();
    let mut wfc = PawWavecar::from_file(format!("{dir}/WAVECAR")).unwrap();

    let ikpt = 1_usize;
    let iband = 10_usize;

    let kvec = wfc.kpoints[ikpt - 1].kvec;
    let bcell = wfc.bcell();
    let encut = wfc.header.encut;

    // Compute projections
    let nonlq = Nonlq::new(&poscar, &pawpot, encut, &kvec, &bcell, false);
    let coeffs = wfc.read_band_coeff(1, ikpt, iband).unwrap();
    let beta = nonlq.proj(&coeffs);

    // Load reference (computed by VaspBandUnfolding paw.py)
    let ref_beta = read_npy_complex128(&format!("{dir}/beta.npy"));

    assert_eq!(beta.len(), ref_beta.len(), "beta length mismatch");

    let max_diff = beta
        .iter()
        .zip(ref_beta.iter())
        .map(|(b, r)| {
            let d = b - r;
            (d.re * d.re + d.im * d.im).sqrt()
        })
        .fold(0.0_f64, f64::max);

    println!("max |beta - ref| = {max_diff:.6e}");

    // The Python reference itself differs from VASP by ~1.4e-5.
    // Our Rust result should match Python to within ~1e-4 (f32 input precision).
    assert!(
        max_diff < 1e-3,
        "max diff {max_diff:.4e} exceeds tolerance; first 5 beta: {:?}",
        &beta[..5.min(beta.len())]
    );
}

#[test]
fn test_nonlq_lmmax() {
    let dir = "tests/pawpot/projectors_lreal_false";
    let poscar = PawPoscar::from_file(format!("{dir}/POSCAR")).unwrap();
    let pawpot = PawPotcar::from_file(&format!("{dir}/POTCAR")).unwrap();
    let wfc = PawWavecar::from_file(format!("{dir}/WAVECAR")).unwrap();

    let kvec = wfc.kpoints[0].kvec;
    let bcell = wfc.bcell();
    let nonlq = Nonlq::new(&poscar, &pawpot, wfc.header.encut, &kvec, &bcell, false);

    // MoS₂: Mo has lmmax=18, S has lmmax=8, total = 18 + 8 + 8 = 34
    assert_eq!(nonlq.total_lmmax, 34, "total lmmax = {}", nonlq.total_lmmax);
    assert_eq!(nonlq.nplw, 2363, "nplw = {}", nonlq.nplw);
}


#[test]
fn test_print_nablaij() {
    let pp = rsgrad::pawpot::PawPotcar::from_file(
        "tests/pawpot/projectors_lreal_false/POTCAR"
    ).unwrap();
    
    for i in 0..pp.len() {
        let pot = pp.get(i).unwrap();
        let nab = pot.get_nablaij();
        let lmmax = pot.lmmax();
        eprintln!("\n=== {} lmmax={} ===", pot.symbol(), lmmax);
        for (alpha, name) in ["Nx", "Ny", "Nz"].iter().enumerate() {
            eprintln!("nablaij[{}] ({}):", alpha, name);
            for row in &nab[alpha] {
                let strs: Vec<String> = row.iter().map(|v| format!("{:9.6}", v)).collect();
                eprintln!("  [{}]", strs.join(", "));
            }
        }
    }
}

#[test]
fn test_print_beta_projections() {
    let dir = "tests/pawpot/projectors_lreal_false";
    let poscar  = rsgrad::pawpot::PawPoscar::from_file(format!("{}/POSCAR", dir)).unwrap();
    let pawpot  = rsgrad::pawpot::PawPotcar::from_file(&format!("{}/POTCAR", dir)).unwrap();
    let mut wav = rsgrad::pawpot::PawWavecar::from_file(format!("{}/WAVECAR", dir)).unwrap();
    
    let encut  = wav.header.encut;
    let kvec   = wav.kpoints[0].kvec;
    let bcell  = wav.bcell();
    
    let nonlq = rsgrad::pawpot::Nonlq::new(&poscar, &pawpot, encut, &kvec, &bcell, false);
    
    for iband in [1usize, 6, 7] {
        let coeff = wav.read_band_coeff(1, 1, iband).unwrap();
        let beta = nonlq.proj(&coeff);
        eprintln!("Band {}: beta[0..8] = {:?}", iband, &beta[..8.min(beta.len())]);
        eprintln!("  S1 part [18..23] = {:?}", &beta[18..23.min(beta.len())]);
    }
}

/// Debug test: print raw dp_ps and PAW correction for bands 6→7 at kpt 1.
/// Compare with VBU reference to diagnose units/formula issues.
#[test]
fn test_debug_tdm_dp_components() {
    use rsgrad::pawpot::gvec::gvectors;

    let dir = "tests/pawpot/projectors_lreal_false";
    let poscar = PawPoscar::from_file(format!("{dir}/POSCAR")).unwrap();
    let pawpot = PawPotcar::from_file(&format!("{dir}/POTCAR")).unwrap();
    let mut wav = PawWavecar::from_file(format!("{dir}/WAVECAR")).unwrap();

    let ikpt  = 1usize;
    let iband = 6usize;
    let jband = 7usize;

    let encut = wav.header.encut;
    let bcell = wav.bcell();
    let kvec  = wav.kpoints[ikpt - 1].kvec;
    let ngrid = wav.header.ngrid;

    let nonlq = Nonlq::new(&poscar, &pawpot, encut, &kvec, &bcell, false);

    let (_, gk, _) = gvectors(&bcell, encut, &kvec, &ngrid);
    let nplw = gk.len();

    let ci = wav.read_band_coeff(1, ikpt, iband).unwrap();
    let cj = wav.read_band_coeff(1, ikpt, jband).unwrap();

    assert_eq!(ci.len(), nplw, "ci len mismatch: {} vs {}", ci.len(), nplw);

    // dp_ps = sum_G conj(Cg_j) * Cg_i * (G+k)_alpha
    let mut dp_ps = [C128::new(0.0, 0.0); 3];
    for ig in 0..nplw {
        let cjc = cj[ig].conj();
        for (alpha, dp_ps_alpha) in dp_ps.iter_mut().enumerate() {
            *dp_ps_alpha += cjc * ci[ig] * gk[ig][alpha];
        }
    }
    eprintln!("dp_ps[x,y,z] = {:.6e} {:.6e} {:.6e}", dp_ps[0], dp_ps[1], dp_ps[2]);
    eprintln!("  |dp_ps| [x,y,z] = {:.6e} {:.6e} {:.6e}",
        dp_ps[0].norm(), dp_ps[1].norm(), dp_ps[2].norm());

    // PAW correction
    let bi = nonlq.proj(&ci);
    let bj = nonlq.proj(&cj);

    let nablaij: Vec<Vec<Vec<Vec<f64>>>> = (0..pawpot.len())
        .map(|it| pawpot.get(it).unwrap().get_nablaij())
        .collect();

    let element_idx = poscar.element_idx();
    let lmmax_per_type: Vec<usize> = (0..pawpot.len())
        .map(|it| pawpot.get(it).unwrap().lmmax())
        .collect();

    let mut corr_total = [C128::new(0.0, 0.0); 3];
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
            corr_total[alpha] += corr;
        }
        nproj += lmmax;
    }
    eprintln!("PAW corr (before -i) [x,y,z] = {:.6e} {:.6e} {:.6e}",
        corr_total[0], corr_total[1], corr_total[2]);
    eprintln!("  |corr| [x,y,z] = {:.6e} {:.6e} {:.6e}",
        corr_total[0].norm(), corr_total[1].norm(), corr_total[2].norm());

    // -i * correction added to dp
    let mut dp_ae = dp_ps;
    for alpha in 0..3 {
        dp_ae[alpha] += C128::new(corr_total[alpha].im, -corr_total[alpha].re);
    }
    eprintln!("dp_ae (after PAW corr) = {:.6e} {:.6e} {:.6e}", dp_ae[0], dp_ae[1], dp_ae[2]);
    eprintln!("  |dp_ae| [x,y,z] = {:.6e} {:.6e} {:.6e}",
        dp_ae[0].norm(), dp_ae[1].norm(), dp_ae[2].norm());

    let ei = wav.kpoints[ikpt - 1].eigenvalues[iband - 1];
    let ej = wav.kpoints[ikpt - 1].eigenvalues[jband - 1];
    let de = ej - ei;
    eprintln!("E6={:.4} E7={:.4} dE={:.4}", ei, ej, de);

    const AUTOA: f64 = 0.529177249;
    const RYTOEV: f64 = 13.605826;
    const AUTDEBYE_VBU: f64 = 2.541746;
    let fac = -1.0 / (de / (2.0 * RYTOEV)) * AUTOA * AUTDEBYE_VBU;
    let dipole_ae: [C128; 3] = std::array::from_fn(|alpha| {
        C128::new(dp_ae[alpha].im * fac, -dp_ae[alpha].re * fac)
    });
    eprintln!("dipole_ae (AUTDEBYE=2.541746) [x,y,z] = {:.6e} {:.6e} {:.6e}",
        dipole_ae[0], dipole_ae[1], dipole_ae[2]);
    eprintln!("  |dipole| [x,y,z] Debye = {:.4} {:.4} {:.4}",
        dipole_ae[0].norm(), dipole_ae[1].norm(), dipole_ae[2].norm());

    let fac_wrong = -1.0 / (de / (2.0 * RYTOEV)) * AUTOA * 4.803242_f64;
    let dipole_ae_wrong: [C128; 3] = std::array::from_fn(|alpha| {
        C128::new(dp_ae[alpha].im * fac_wrong, -dp_ae[alpha].re * fac_wrong)
    });
    eprintln!("dipole_ae (AUTDEBYE=4.803242) [x,y,z] = {:.6e} {:.6e} {:.6e}",
        dipole_ae_wrong[0], dipole_ae_wrong[1], dipole_ae_wrong[2]);
    eprintln!("  |dipole| [x,y,z] Debye = {:.4} {:.4} {:.4}",
        dipole_ae_wrong[0].norm(), dipole_ae_wrong[1].norm(), dipole_ae_wrong[2].norm());

    // Per-row breakdown for Mo (atom 0, lmmax=18)
    eprintln!("\n--- Per-atom correction debug ---");
    let mut nproj2 = 0usize;
    for (iatom, &ntype) in element_idx.iter().enumerate() {
        let lmmax = lmmax_per_type[ntype];
        let nab = &nablaij[ntype];
        eprintln!("Atom {} (type={}, lmmax={}, nproj={}):", iatom, ntype, lmmax, nproj2);
        for (alpha, nab_alpha) in nab.iter().enumerate().take(1) { // just x
            let mut corr = C128::new(0.0, 0.0);
            for m in 0..lmmax {
                let mut tmp = C128::new(0.0, 0.0);
                for n in 0..lmmax {
                    tmp += nab_alpha[m][n] * bi[nproj2 + n];
                }
                let contrib = bj[nproj2 + m].conj() * tmp;
                if contrib.norm() > 1e-8 {
                    eprintln!("  m={}: bj[{}]={:.4e} tmp={:.4e} contrib={:.4e}",
                        m, nproj2+m, bj[nproj2+m], tmp, contrib);
                }
                corr += contrib;
            }
            eprintln!("  corr[alpha={}]={:.6e}", alpha, corr);
        }
        // Print bi and bj for this atom
        eprintln!("  bi[{}..{}] norms = {:?}", nproj2, nproj2+lmmax,
            bi[nproj2..nproj2+lmmax].iter().map(|v| format!("{:.3e}", v.norm())).collect::<Vec<_>>());
        eprintln!("  bj[{}..{}] norms = {:?}", nproj2, nproj2+lmmax,
            bj[nproj2..nproj2+lmmax].iter().map(|v| format!("{:.3e}", v.norm())).collect::<Vec<_>>());
        nproj2 += lmmax;
    }
}

/// Debug test: print sum of |sg_ae|², |sg_ps|², and |sg_ae - sg_ps|² for band 1, kpt 1.
/// Uses MoS2 test data (WAVECAR, POSCAR, POTCAR in projectors_lreal_false/).
#[test]
fn test_debug_sg_ae_ps_norm() {
    let dir = "tests/pawpot/projectors_lreal_false";
    let poscar = PawPoscar::from_file(format!("{dir}/POSCAR")).unwrap();
    let pawpot = PawPotcar::from_file(&format!("{dir}/POTCAR")).unwrap();
    let wavecar = PawWavecar::from_file(format!("{dir}/WAVECAR")).unwrap();

    let ikpt = 1usize;
    let iband = 1usize;
    let ispin = 1usize;

    let mut ae = VaspAeWfc::new(wavecar, &poscar, &pawpot, ikpt, -4.0).unwrap();

    let (sg_ae, sg_ps) = ae.compute_sg_ae_ps(ispin, iband).unwrap();

    let n_ae = sg_ae.len();

    let sum_ae_sq: f64 = sg_ae.iter().map(|v| v.norm_sqr()).sum();
    let sum_ps_sq: f64 = sg_ps.iter().map(|v| v.norm_sqr()).sum();
    let sum_diff_sq: f64 = sg_ae.iter()
        .zip(sg_ps.iter())
        .map(|(ae, ps)| (ae - ps).norm_sqr())
        .sum();

    println!("MoS2 band={iband} kpt={ikpt}: n_ae_gvecs={n_ae}");
    println!("  sum |sg_ae[ig]|^2      = {sum_ae_sq:.8e}");
    println!("  sum |sg_ps[ig]|^2      = {sum_ps_sq:.8e}");
    println!("  sum |sg_ae-sg_ps|^2    = {sum_diff_sq:.8e}");
}
