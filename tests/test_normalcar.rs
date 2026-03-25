// Integration test: write NormalCAR and compare cproj against reference (cproj.npy).

use std::convert::TryInto;

use rsgrad::pawpot::{write_normalcar, PawPotcar, PawPoscar, PawWavecar};
use num::complex::Complex;

type C128 = Complex<f64>;

/// Minimal .npy reader for complex128 arrays.
fn read_npy_complex128(path: &str) -> Vec<C128> {
    let bytes = std::fs::read(path).expect("read npy file");
    assert_eq!(&bytes[0..6], b"\x93NUMPY", "bad npy magic");
    let major = bytes[6];
    let header_len = if major == 1 {
        u16::from_le_bytes([bytes[8], bytes[9]]) as usize
    } else {
        u32::from_le_bytes([bytes[8], bytes[9], bytes[10], bytes[11]]) as usize
    };
    let header_offset = if major == 1 { 10 } else { 12 };
    let data_start = header_offset + header_len;
    let raw = &bytes[data_start..];
    let n = raw.len() / 16;
    (0..n)
        .map(|i| {
            let re = f64::from_le_bytes(raw[16 * i..16 * i + 8].try_into().unwrap());
            let im = f64::from_le_bytes(raw[16 * i + 8..16 * i + 16].try_into().unwrap());
            C128::new(re, im)
        })
        .collect()
}

/// Read all band projectors from a NormalCAR file, skipping the 4 header records.
fn read_normalcar_cproj(path: &str) -> Vec<Vec<C128>> {
    let bytes = std::fs::read(path).expect("read NormalCAR");
    let mut offset = 0usize;

    let read_rec = |bytes: &[u8], off: &mut usize| -> Vec<u8> {
        let marker = i32::from_le_bytes(bytes[*off..*off + 4].try_into().unwrap()) as usize;
        let payload = bytes[*off + 4..*off + 4 + marker].to_vec();
        let end = i32::from_le_bytes(bytes[*off + 4 + marker..*off + 8 + marker].try_into().unwrap()) as usize;
        assert_eq!(marker, end, "record marker mismatch at offset {}", *off);
        *off += 8 + marker;
        payload
    };

    // Skip 4 header records: rec1(lmdim/nions/nrspinors), rec2(cqij), rec3(nprod/npro/ntyp), type0
    for _ in 0..4 {
        read_rec(&bytes, &mut offset);
    }

    // Read band records
    let mut records = Vec::new();
    while offset < bytes.len() {
        let payload = read_rec(&bytes, &mut offset);
        let n = payload.len() / 16;
        let cproj: Vec<C128> = (0..n)
            .map(|i| {
                let re = f64::from_le_bytes(payload[16 * i..16 * i + 8].try_into().unwrap());
                let im = f64::from_le_bytes(payload[16 * i + 8..16 * i + 16].try_into().unwrap());
                C128::new(re, im)
            })
            .collect();
        records.push(cproj);
    }
    records
}

#[test]
fn test_write_normalcar_lreal_false() {
    let dir = "tests/pawpot/projectors_lreal_false";
    let out = "/tmp/test_NormalCAR_lreal_false";

    let poscar = PawPoscar::from_file(format!("{dir}/POSCAR")).unwrap();
    let pawpot = PawPotcar::from_file(&format!("{dir}/POTCAR")).unwrap();
    let mut wavecar = PawWavecar::from_file(format!("{dir}/WAVECAR")).unwrap();

    write_normalcar(out, &poscar, &pawpot, &mut wavecar).unwrap();

    // Load reference cproj (flat: nkpts*nbands × npro_tot)
    let ref_flat = read_npy_complex128(&format!("{dir}/cproj.npy"));
    let npro_tot = ref_flat.len() / (wavecar.header.nkpts * wavecar.header.nbands);

    // Read our output
    let our_records = read_normalcar_cproj(out);
    let n_records = wavecar.header.nspin * wavecar.header.nkpts * wavecar.header.nbands;
    assert_eq!(our_records.len(), n_records,
        "expected {} band records, got {}", n_records, our_records.len());
    assert_eq!(our_records[0].len(), npro_tot,
        "expected {} projectors per band, got {}", npro_tot, our_records[0].len());

    // Compare against reference (cproj.npy is flat: shape (nkpts*nbands, npro_tot))
    let max_diff = our_records
        .iter()
        .flatten()
        .zip(ref_flat.iter())
        .map(|(a, b)| (a - b).norm())
        .fold(0.0_f64, f64::max);

    println!("max |cproj - ref| = {:.4e}", max_diff);
    assert!(
        max_diff < 1e-3,
        "max diff {:.4e} exceeds tolerance 1e-3",
        max_diff
    );
}
