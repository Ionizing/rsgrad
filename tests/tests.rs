use std::path::PathBuf;
use std::io;
use rsgrad::outcar::Outcar;

// #[macro_export]
macro_rules! get_fpath_in_current_dir {
    ($fname:expr) => {{
        let mut path = PathBuf::from(file!());
        path.pop();
        path.push($fname);
        path
    }}
}

#[test]
fn test_normal_outcar() -> io::Result<()> {
    let fname = get_fpath_in_current_dir!("OUTCAR_multiple_ionic_steps");
    let outcar = Outcar::from_file(&fname)?;

    assert_eq!(outcar.lsorbit, false);
    assert_eq!(outcar.ispin, 1);
    assert_eq!(outcar.ibrion, 1);
    assert_eq!(outcar.nions, 32);
    assert_eq!(outcar.nkpts, 20);
    assert_eq!(outcar.nbands, 81);
    assert_eq!(outcar.efermi, 2.9331);
    assert_eq!(outcar.cell, [[7.519999981,         0.0,         0.0],
                             [        0.0, 7.519999981,         0.0],
                             [        0.0,         0.0, 7.519999981]]);
    assert_eq!(outcar.ext_pressure, vec![-18.05, 21.53, -2.72, -5.24, -0.30]);
    assert_eq!(outcar.ions_per_type, vec![32]);
    assert_eq!(outcar.ion_types, vec!["C"]);
    assert_eq!(outcar.ion_masses, vec![12.011; 32]);
    assert_eq!(outcar.ion_iters.len(), 5);
    assert_eq!(outcar.vib, None);
    outcar.ion_iters.iter()
                    .zip(vec![14i32, 8, 7, 8, 7].iter())
                    .for_each(|(x, y)| assert_eq!(&x.nscf, y));

    outcar.ion_iters.iter()
                    .zip(vec![-253.61858820,
                              -253.61023247,
                              -253.61629491,
                              -253.58960211,
                              -253.64363797].iter())
                    .for_each(|(x, y)| assert_eq!(&x.toten, y));

    outcar.ion_iters.iter()
                    .zip(vec![-253.61858820,
                              -253.61023247,
                              -253.61629491,
                              -253.58960211,
                              -253.64363797].iter())
                    .for_each(|(x, y)| assert_eq!(&x.toten, y));
    Ok(())
}


#[test]
fn test_ispin2_outcar() -> io::Result<()> {
    let fname = get_fpath_in_current_dir!("OUTCAR_ispin2");
    let outcar = Outcar::from_file(&fname)?;
    Ok(())
}


#[test]
fn test_ncl_outcar() -> io::Result<()> {
    let fname = get_fpath_in_current_dir!("OUTCAR_ncl");
    let outcar = Outcar::from_file(&fname)?;
    Ok(())
}
