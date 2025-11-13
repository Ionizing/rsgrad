//use std::io::Read;
use std::io::{
    Seek,
    SeekFrom,
};
use std::fs;
use std::path::Path;

use anyhow::{self, Result, Context};
use ndarray as na;
use byteorder::{
    LittleEndian,
    ReadBytesExt,
};
use ndrustfft::Complex;

#[allow(non_camel_case_types)]
type c64 = Complex<f64>;
#[allow(non_camel_case_types)]
type c32 = Complex<f32>;


/// Read NormalCAR data, ikpoint count from 1, then returns the projector coefficients CPROJ
/// shape(CPROJ) = (2, nbands, nproj)
pub fn read_normalcar<P>(fname: P, nbands: usize, nkpoints: usize, ikpoint: usize) -> Result<(na::Array3<c64> /* cproj */, usize /* nproj */)>
where P: AsRef<Path> {
    let mut f = fs::File::open(&fname).context(format!("Failed to open file {:?}.", fname.as_ref()))?;

    // rec_l, lmdim, nions, nrspinors, rec_r
    let mut buf = [0i32; 5];
    f.read_i32_into::<LittleEndian>(&mut buf).unwrap();
    let [rec_l, lmdim, nions, nrspinors, rec_r] = buf;
    assert_eq!(rec_l, rec_r, "Invalid record length.");

    // skip cqij, read() rec_l, cqij(lmdim, lmdim, nions, nrspinors), rec_r
    let rec_l = f.read_i32::<LittleEndian>().unwrap();
    f.seek(SeekFrom::Current((8 * lmdim * lmdim * nions * nrspinors) as i64)).unwrap();
    let rec_r = f.read_i32::<LittleEndian>().unwrap();
    assert_eq!(rec_l, rec_r, "Invalid record length.");


    // rec_l, nprod, npro, ntyp, rec_r
    let mut buf = [0i32; 5];
    f.read_i32_into::<LittleEndian>(&mut buf).unwrap();
    let [rec_l, _nprod, _npro, ntyp, rec_r] = buf;
    assert_eq!(rec_l, rec_r, "Invalid record length.");

    // skip [rec_l, lmmax, nityp, rec_r] * ntyp
    f.seek(SeekFrom::Current(4 * (ntyp * 4) as i64)).unwrap();


    let rec_l = f.read_i32::<LittleEndian>().unwrap();
    let nproj = rec_l / 16;
    f.seek(SeekFrom::Current(-4)).unwrap();
    let mut cproj = na::Array3::<c64>::zeros((2, nbands as usize, nproj as usize));
    let mut buf = vec![0.0f64; nproj as usize * 2];
    for ispin in 0 .. 2 {
        for ikpt in 0 .. nkpoints {
            for iband in 0 .. nbands {
                let rec_l = f.read_i32::<LittleEndian>().unwrap();
                
                if ikpt + 1 == ikpoint {
                    f.read_f64_into::<LittleEndian>(&mut buf).unwrap();
                    for iproj in 0 .. nproj {
                        cproj[(ispin as usize, iband as usize, iproj as usize)] = c64::new(buf[2 * iproj as usize], buf[2 * iproj as usize + 1]);
                    }
                } else {
                    f.seek(SeekFrom::Current(16 * nproj as i64)).unwrap();
                }

                let rec_r = f.read_i32::<LittleEndian>().unwrap();
                assert_eq!(rec_l, rec_r, "Invalid record length.");
            }
        }
    }

    Ok((cproj, nproj as usize))
}


/// Read SocCar, returns the SOC matrix.
/// shape(SOC) = (4, nproj, nproj)
/// where: SOC[0, .., ..] = up to up
///        SCO[1, .., ..] = up to dn
///        SCO[2, .., ..] = dn to up
///        SCO[3, .., ..] = dn to dn
pub fn read_soccar<P>(fname: P, nproj: usize) -> Result<na::Array3<c64> /* soc */>
where P: AsRef<Path> {
    let txt = fs::read_to_string(&fname).context(format!("Failed to open file: {:?}.", fname.as_ref()))?;

    let v = txt.split_ascii_whitespace()
        .map(|x| {
            x.parse::<f64>()
                .with_context(|| format!("Cannot parse {} as f64.", x))
                .unwrap()
        })
        .collect::<Vec<f64>>();

    //let nproj = ((v.len() / 8) as f64).sqrt().round() as usize;
    anyhow::ensure!(nproj * nproj * 4 * 2 == v.len(), "Invalid SocCar length.");

    let ret = na::Array1::from_vec(vec_to_complex(v));
    Ok(ret.into_shape((4, nproj, nproj))?)
}


/// Calculate band-to-band spin-orbit coupling, needs reading SocCar and NormalCar
///
/// ikpoint counts from 1
///
/// Hmm layout: [uu, ud, du, dd]
#[allow(non_snake_case)]
pub fn calc_hmm<P>(runpath: P, nbands: usize, nkpoints: usize, ikpoint: usize) -> Result<na::Array3<c64> /* Hmm */>
where P: AsRef<Path> {
    let normalcar_fname = runpath.as_ref().join("NormalCAR");
    let soccar_fname = runpath.as_ref().join("SocCar");

    let (cproj, nproj) = read_normalcar(normalcar_fname, nbands, nkpoints, ikpoint)?;
    let soccar = read_soccar(soccar_fname, nproj)?;

    Ok(calc_hmm_helper(&cproj, &soccar))
}


pub fn calc_hmm_helper(cproj: &na::Array3<c64>, soccar: &na::Array3<c64>) -> na::Array3<c64> {
    let cproj_shape = cproj.shape();
    assert_eq!(cproj_shape[0], 2);
    let nbands = cproj_shape[1];
    let nproj  = cproj_shape[2];

    let soccar_shape = soccar.shape();
    assert_eq!(4, soccar_shape[0]);
    assert_eq!(nproj, soccar_shape[1]);
    assert_eq!(nproj, soccar_shape[2]);

    let mut hmm = na::Array3::<c64>::zeros((4, nbands, nbands));

    // CPROJ[0, .., ..]^H * SOCCAR * CPROJ[0, .., ..]
    hmm.index_axis_mut(na::Axis(0), 0).assign(
        &cproj.index_axis(na::Axis(0), 0).mapv(|v| v.conj())
            .dot(&soccar.index_axis(na::Axis(0), 0))
            .dot(&cproj.index_axis(na::Axis(0), 0).t())
    );

    hmm.index_axis_mut(na::Axis(0), 3).assign(
        &cproj.index_axis(na::Axis(0), 1).mapv(|v| v.conj())
            .dot(&soccar.index_axis(na::Axis(0), 3))
            .dot(&cproj.index_axis(na::Axis(0), 1).t())
    );

    hmm.index_axis_mut(na::Axis(0), 2).assign(
        &cproj.index_axis(na::Axis(0), 1).mapv(|v| v.conj())
            .dot(&soccar.index_axis(na::Axis(0), 2))
            .dot(&cproj.index_axis(na::Axis(0), 0).t())
    );

    hmm.index_axis_mut(na::Axis(0), 1).assign(
        &cproj.index_axis(na::Axis(0), 0).mapv(|v| v.conj())
            .dot(&soccar.index_axis(na::Axis(0), 1))
            .dot(&cproj.index_axis(na::Axis(0), 1).t())
    );

    hmm
}


// https://stackoverflow.com/a/54188098/8977923
fn vec_to_complex(mut buffer: Vec<f64>) -> Vec<c64> {
    unsafe {
        buffer.shrink_to_fit();
        
        let ptr = buffer.as_mut_ptr() as *mut c64;
        let len = buffer.len();
        let cap = buffer.capacity();

        assert!(len % 2 == 0);
        assert!(cap % 2 == 0);

        std::mem::forget(buffer);

        Vec::from_raw_parts(ptr, len / 2, cap / 2)
    }
}


// https://stackoverflow.com/a/54188098/8977923
fn slice_to_complex(buffer: &[f64]) -> &[c64] {
    unsafe {
        let ptr = buffer.as_ptr() as *mut c64;
        let len = buffer.len();

        assert!(len % 2 == 0);
        std::slice::from_raw_parts(ptr, len / 2)
    }
}


#[cfg(test)]
mod tests {
    use super::*;

    #[test]
    fn test_vec_to_complex() {
        let v = vec![1.0, 2.0, 5.0, 6.0];
        let cv = vec_to_complex(v);
        assert_eq!(cv, vec![ c64::new(1.0, 2.0), c64::new(5.0, 6.0) ]);

        let n = 65536usize;
        let v = (0 .. n).map(|x| x as f64).collect::<Vec<f64>>();
        let cv = vec_to_complex(v);
        let cv_expected = (0 .. n/2).map(|x| c64::new(x as f64 * 2.0, x as f64 * 2.0 + 1.0)).collect::<Vec<c64>>();
        assert_eq!(cv, cv_expected);
    }


    #[test]
    fn test_slice_to_complex() {
        let n = 65536usize;
        let v = (0 .. n).map(|x| x as f64).collect::<Vec<f64>>();
        let cv = slice_to_complex(&v);
        let cv_expected = (0 .. n/2).map(|x| c64::new(x as f64 * 2.0, x as f64 * 2.0 + 1.0)).collect::<Vec<c64>>();
        assert_eq!(cv, cv_expected);
    }

    #[test]
    fn test_read_normalcar() {
        let nbands = 208usize;
        let nkpoints = 14usize;
        let ikpoint = 1usize;
        let (cproj, nproj) = read_normalcar("tests/NormalCAR", nbands, nkpoints, ikpoint).unwrap();
        assert_eq!(nproj, 576);
        assert_eq!(cproj[(0, 0, 0)], c64::new(-0.0019005957560536114, 0.0043363155065891225));
        assert_eq!(cproj[(0, 0, 1)], c64::new(-0.00033174315808717357, 0.0007570752911495605));

        assert_eq!(cproj[(0, 1, 0)], c64::new(-0.0008910190010085349, 0.0008318299869942524));
        assert_eq!(cproj[(0, 1, 1)], c64::new(-0.0006250675689409707, 0.0005841185595735311));
    }

    #[test]
    fn test_read_soccar() {
        //let nbands = 208usize;
        //let nkpoints = 14usize;
        //let ikpoint = 1usize;
        let nproj = 576;
        let soccar = read_soccar("tests/SocCar", 576).unwrap();
        assert_eq!(soccar[(3, nproj-3, nproj-1)], c64::new(0.0, -0.1501573));
        assert_eq!(soccar[(3, nproj-1, nproj-3)], c64::new(0.0,  0.1501573));
    }

    #[test]
    fn test_calc_hmm() {
        let nbands = 208usize;
        let nkpoints = 14usize;
        let ikpoint = 1usize;
        let hmm = calc_hmm("tests", nbands, nkpoints, ikpoint).unwrap();
        assert_eq!(hmm[(0, 0, 0)], c64::new(-5.5624424817112524e-12, -7.940933880509066e-22));
        assert_eq!(hmm[(1, 207, 0)], c64::new(0.0002937153288168673, 3.179271745971495e-5));
        assert_eq!(hmm[(0, 207, 207)], c64::new(-3.608022704651101e-5, 2.1277467635028025e-18));
    }
}
