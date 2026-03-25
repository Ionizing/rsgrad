/// Benchmark: Compare pseudo-TDM vs AE-TDM from rsgrad
use rsgrad::pawpot::*;
use rsgrad::pawpot::gvec::gvectors;
use num::complex::Complex;

type C128 = Complex<f64>;

#[test]
fn bench_tdm_accuracy() {
    let dir = "tests/pawpot/projectors_lreal_false";

    println!("\n{}", "=".repeat(100));
    println!("TDM Accuracy Benchmark: Pseudo vs AE-TDM");
    println!("{}\n", "=".repeat(100));

    let poscar = PawPoscar::from_file(format!("{}/POSCAR", dir)).unwrap();
    let pawpot = PawPotcar::from_file(&format!("{}/POTCAR", dir)).unwrap();
    let mut wavecar = PawWavecar::from_file(format!("{}/WAVECAR", dir)).unwrap();

    let encut = wavecar.header.encut;
    let bcell = wavecar.bcell();
    let nbands = wavecar.header.nbands;
    let nkpts = wavecar.header.nkpts;

    println!("System:");
    println!("  File: {}/WAVECAR", dir);
    println!("  nbands: {}, nkpts: {}, encut: {:.2} eV\n", nbands, nkpts, encut);

    let mut tdm_stats = TdmStats::new();

    // First k-point only (to save time)
    let ikpt = 1;
    let kvec = wavecar.kpoints[ikpt - 1].kvec;

    println!("Computing TDM for k-point 1 (kvec = [{:.4}, {:.4}, {:.4}])...\n",
             kvec[0], kvec[1], kvec[2]);

    // Compute AE-TDM using rsgrad
    let ibands: Vec<usize> = (1..=std::cmp::min(5, nbands)).collect();
    let jbands: Vec<usize> = (2..=std::cmp::min(8, nbands)).collect();

    let entries = compute_tdm_paw(
        &poscar,
        &pawpot,
        &mut wavecar,
        1,       // ispin
        ikpt,    // ikpt
        &ibands,
        &jbands,
        false,   // not lpseudo
    ).unwrap_or_default();

    println!("{:<8} {:<8} {:<15} {:<15} {:<15} {:<12}",
             "i", "j", "|P_pseudo|", "|P_AE|", "|ΔP|", "PAW_corr%");
    println!("{}", "-".repeat(80));

    // For comparison, compute pseudo-TDM manually
    let volume = poscar.volume().abs();
    let inv_sqrt_vol = 1.0 / volume.sqrt();

    let ngrid = wavecar.header.ngrid;
    let (gvec_all, _, _) = gvectors(&bcell, encut, &kvec, &ngrid);

    for iband in 1..=std::cmp::min(5, nbands) {
        let coeffs_i = wavecar.read_band_coeff(1, ikpt, iband).unwrap();

        for jband in (iband + 1)..=std::cmp::min(iband + 3, nbands) {
            let coeffs_j = wavecar.read_band_coeff(1, ikpt, jband).unwrap();

            // Compute pseudo-TDM (reciprocal space)
            let mut p_ps = C128::new(0.0, 0.0);
            let gvec = &gvec_all;

            for (ig, &g) in gvec.iter().enumerate() {
                if ig >= coeffs_i.len() || ig >= coeffs_j.len() {
                    break;
                }

                let gk = [
                    kvec[0] + g[0] as f64,
                    kvec[1] + g[1] as f64,
                    kvec[2] + g[2] as f64,
                ];

                let gk_cart = [
                    gk[0] * bcell[0][0] + gk[1] * bcell[1][0] + gk[2] * bcell[2][0],
                    gk[0] * bcell[0][1] + gk[1] * bcell[1][1] + gk[2] * bcell[2][1],
                    gk[0] * bcell[0][2] + gk[1] * bcell[1][2] + gk[2] * bcell[2][2],
                ];

                let conj_i = coeffs_i[ig].conj();
                let gk_vec = C128::new(0.0, 1.0) * C128::new(
                    gk_cart[0] + gk_cart[1] + gk_cart[2],
                    0.0,
                );

                p_ps += conj_i * coeffs_j[ig] * gk_vec;
            }

            let p_pseudo = p_ps * inv_sqrt_vol;
            let p_pseudo_norm = p_pseudo.norm();

            // Find corresponding AE-TDM
            if let Some(entry) = entries.iter().find(|e| e.iband == iband && e.jband == jband) {
                let p_ae_norm = entry.dipole[0].norm();
                let delta = (p_ae_norm - p_pseudo_norm).abs();
                let paw_corr_pct = if p_pseudo_norm > 1e-10 {
                    (delta / p_pseudo_norm) * 100.0
                } else {
                    0.0
                };

                println!("{:<8} {:<8} {:<15.8} {:<15.8} {:<15.8} {:<11.2}%",
                         iband, jband, p_pseudo_norm, p_ae_norm, delta, paw_corr_pct);

                tdm_stats.add(p_pseudo_norm, p_ae_norm, delta);
            } else {
                println!("{:<8} {:<8} {:<15.8} {:15} {:15} {:12}",
                         iband, jband, p_pseudo_norm, "(no AE)", "(no AE)", "(no AE)");
            }
        }
    }

    println!("\n{}", "-".repeat(80));
    tdm_stats.print_summary();

    // Verify reasonable values
    if tdm_stats.count > 0 {
        let avg_pseudo = tdm_stats.total_pseudo / (tdm_stats.count as f64);
        let avg_ae = tdm_stats.total_ae / (tdm_stats.count as f64);

        println!("\nVerification:");
        println!("  Pseudo-TDM values typically 0.1-0.5 (atomic units)");
        println!("  PAW correction typically 1-5%");
        println!("  Avg |P_pseudo| = {:.6}", avg_pseudo);
        println!("  Avg |P_AE| = {:.6}", avg_ae);

        assert!(avg_pseudo > 1e-6, "Pseudo-TDM seems too small: {:.2e}", avg_pseudo);
    }

    println!("{}\n", "=".repeat(100));
}

struct TdmStats {
    count: usize,
    total_pseudo: f64,
    total_ae: f64,
    total_delta: f64,
}

impl TdmStats {
    fn new() -> Self {
        TdmStats {
            count: 0,
            total_pseudo: 0.0,
            total_ae: 0.0,
            total_delta: 0.0,
        }
    }

    fn add(&mut self, p_pseudo: f64, p_ae: f64, delta: f64) {
        self.count += 1;
        self.total_pseudo += p_pseudo;
        self.total_ae += p_ae;
        self.total_delta += delta;
    }

    fn print_summary(&self) {
        if self.count == 0 {
            println!("No transitions computed");
            return;
        }

        let avg_pseudo = self.total_pseudo / (self.count as f64);
        let avg_ae = self.total_ae / (self.count as f64);
        let avg_delta = self.total_delta / (self.count as f64);
        let avg_paw_pct = if avg_pseudo > 1e-10 {
            (avg_delta / avg_pseudo) * 100.0
        } else {
            0.0
        };

        println!("Summary ({} transitions):", self.count);
        println!("  Avg |P_pseudo|:      {:.8}", avg_pseudo);
        println!("  Avg |P_AE|:          {:.8}", avg_ae);
        println!("  Avg |ΔP|:            {:.8}", avg_delta);
        println!("  Avg PAW correction:  {:.2}%", avg_paw_pct);
    }
}
