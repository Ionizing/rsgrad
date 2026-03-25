// Integration test for VaspAeWfc.
// Verifies that the AE norm ≈ 1.0 for all bands of CO₂ (gamma-only WAVECAR).

use rsgrad::pawpot::{VaspAeWfc, PawPotcar, PawPoscar, PawWavecar};

#[test]
fn test_ae_norm_co2() {
    let dir = "tests/pawpot/aewfc_co2";

    let poscar = PawPoscar::from_file(format!("{dir}/POSCAR")).unwrap();
    let pawpot = PawPotcar::from_file(&format!("{dir}/POTCAR")).unwrap();
    let wavecar = PawWavecar::from_file(format!("{dir}/WAVECAR")).unwrap();

    let nbands = wavecar.header.nbands;

    let mut ae = VaspAeWfc::new(wavecar, &poscar, &pawpot, 1, -4.0).unwrap();

    assert!(ae.lgam, "Expected gamma-only WAVECAR for CO2");

    let mut max_dev = 0.0f64;
    for iband in 1..=nbands {
        let norm = ae.get_ae_norm(1, iband).unwrap();
        let dev = (norm - 1.0).abs();
        println!("band {:3}: ae_norm = {:.8}", iband, norm);
        if dev > max_dev {
            max_dev = dev;
        }
    }

    println!("max |ae_norm - 1| = {:.2e}", max_dev);
    assert!(
        max_dev < 1e-3,
        "max deviation from 1.0 is {:.4e}, expected < 1e-3",
        max_dev
    );
}

#[test]
fn test_ae_norm_debug_beta() {
    let dir = "tests/pawpot/aewfc_co2";
    let poscar = PawPoscar::from_file(format!("{dir}/POSCAR")).unwrap();
    let pawpot = PawPotcar::from_file(&format!("{dir}/POTCAR")).unwrap();
    let wavecar = PawWavecar::from_file(format!("{dir}/WAVECAR")).unwrap();

    let bcell = wavecar.bcell();
    let kvec = wavecar.kpoints[0].kvec;
    let encut = wavecar.header.encut;

    // Build nonlq with lgam=true
    let nonlq = rsgrad::pawpot::Nonlq::new(&poscar, &pawpot, encut, &kvec, &bcell, true);
    println!("nplw = {}", nonlq.nplw);
    println!("total_lmmax = {}", nonlq.total_lmmax);

    let mut wf2 = PawWavecar::from_file(format!("{dir}/WAVECAR")).unwrap();
    let coeffs = wf2.read_band_coeff(1, 1, 1).unwrap();

    let ps_norm: f64 = coeffs.iter().map(|c| c.norm_sqr()).sum();
    println!("ps_norm = {:.8}", ps_norm);

    let beta = nonlq.proj(&coeffs);
    println!("beta[0:5]: {:?}", &beta[..5]);
}

#[test]
fn test_reciprojs_co2() {
    use rsgrad::pawpot::PawPotcar;
    let pawpot = PawPotcar::from_file(&"tests/pawpot/aewfc_co2/POTCAR").unwrap();
    let pp = pawpot.get(0).unwrap(); // C
    println!("C gmax={}, lmax={}", pp.gmax(), pp.lmax());
    println!("C proj_qgrid[0:5]: {:?}", &pp.proj_qgrid()[..5]);
    println!("C reciprojs[0][0:5]: {:?}", &pp.reciprojs()[0][..5]);
    println!("C ls: {:?}", pp.ls());
    // Check spl_qproj eval at q=0
    let spls = pp.spl_qproj();
    println!("C spl_qproj[0].eval(0) = {:?}", spls[0].eval(0.0));
    println!("C spl_qproj[0].eval(0.6283) = {:?}", spls[0].eval(0.6283));
}

#[test]
fn test_crexp_debug() {
    use rsgrad::pawpot::{PawPoscar, PawPotcar, PawWavecar, Nonlq};
    let dir = "tests/pawpot/aewfc_co2";
    let poscar = PawPoscar::from_file(format!("{dir}/POSCAR")).unwrap();
    println!("positions: {:?}", &poscar.positions[..3]);

    let pawpot = PawPotcar::from_file(&format!("{dir}/POTCAR")).unwrap();
    let wavecar = PawWavecar::from_file(format!("{dir}/WAVECAR")).unwrap();
    let bcell = wavecar.bcell();
    let kvec = wavecar.kpoints[0].kvec;
    let encut = wavecar.header.encut;

    let nonlq = Nonlq::new(&poscar, &pawpot, encut, &kvec, &bcell, true);

    // Print first few G-vectors
    println!("gvec[0:5]: {:?}", &nonlq.gvec[..5]);
    println!("glen[0:5]: {:?}", &nonlq.glen[..5]);

    // Print first few crexp (via proj computation)
    let mut wf2 = PawWavecar::from_file(format!("{dir}/WAVECAR")).unwrap();
    let coeffs = wf2.read_band_coeff(1, 1, 1).unwrap();
    println!("coeffs[0:3]: {:?}", &coeffs[..3]);

    // Manually compute first few qproj contributions
    let pp = pawpot.get(0).unwrap();
    let spls = pp.spl_qproj();
    let gk0 = &nonlq.gk[..5];

    use rsgrad::pawpot::sph_r;
    let ylm0 = sph_r(gk0, 0);
    println!("ylm[l=0][ig=0..4]: {:?}", &ylm0[..4]);

    // Check qproj values
    for ig in 0..5 {
        let q = spls[0].eval(nonlq.glen[ig]).unwrap_or(0.0);
        let y = ylm0[ig][0];
        let vol = poscar.volume().abs();
        let base = q * y / vol.sqrt();
        let factor = if ig == 0 { 1.0 } else { std::f64::consts::SQRT_2 };
        println!("ig={ig}: glen={:.4}, spl={:.4}, ylm={:.4}, qproj_scaled={:.6}",
            nonlq.glen[ig], q, y, base * factor);
    }
}

#[test]
fn test_qij_co2() {
    use rsgrad::pawpot::PawPotcar;
    let pawpot = PawPotcar::from_file(&"tests/pawpot/aewfc_co2/POTCAR").unwrap();
    let pp_c = pawpot.get(0).unwrap();
    let qij_c = pp_c.get_qij();
    println!("C lmmax={}", pp_c.lmmax());
    println!("C ls={:?}", pp_c.ls());
    println!("Q_C diagonal:");
    for i in 0..pp_c.lmmax() {
        print!("  [{i}] = {:.8}", qij_c[i][i]);
    }
    println!();
    println!("Q_C[0][1]={:.8}, Q_C[1][0]={:.8}", qij_c[0][1], qij_c[1][0]);

    let pp_o = pawpot.get(1).unwrap();
    let qij_o = pp_o.get_qij();
    println!("O lmmax={}", pp_o.lmmax());
    println!("Q_O diagonal:");
    for i in 0..pp_o.lmmax() {
        print!("  [{i}] = {:.8}", qij_o[i][i]);
    }
    println!();

    // Compute beta†Qβ for band 1 (C atom only) with known beta from Python
    let beta_c = [2.31227627, -1.23964067, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0_f64];
    let mut qterm_c = 0.0f64;
    for i in 0..8 {
        for j in 0..8 {
            qterm_c += beta_c[i] * qij_c[i][j] * beta_c[j];
        }
    }
    println!("beta_c†Q_C beta_c = {:.8}", qterm_c);
}
