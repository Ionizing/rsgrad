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
                    .for_each(|(x, y)| assert_eq!(&x.toten_z, y));

    assert_eq!(&outcar.ion_iters.last().unwrap().cell, &[[7.494265554, 0.000000000, -0.000000000],
                                                         [0.000000000, 7.494265554, -0.000000000],
                                                         [0.000000000, 0.000000000,  7.494265554]]);

    assert_eq!(outcar.ion_iters.last().unwrap()
               .positions.last().unwrap(), &[5.09135, 5.09135, 1.34421]);
    assert_eq!(outcar.ion_iters.last().unwrap()
               .forces.last().unwrap(), &[-0.000716, -0.000716, -0.000716]);

    assert!(outcar.ion_iters.iter().all(|i| i.magmom.is_none()));
    Ok(())
}


#[test]
fn test_ispin2_outcar() -> io::Result<()> {
    let fname = get_fpath_in_current_dir!("OUTCAR_ispin2");
    let outcar = Outcar::from_file(&fname)?;

    assert_eq!(outcar.lsorbit, false);
    assert_eq!(outcar.ispin, 2);
    assert_eq!(outcar.ibrion, 1);
    assert_eq!(outcar.nions, 3);
    assert_eq!(outcar.nkpts, 41);
    assert_eq!(outcar.nbands, 16);
    assert_eq!(outcar.efermi, -2.2687);
    assert_eq!(outcar.cell, [[ 2.864537506, -1.653841500,  0.000000000],
                             [ 0.000000000,  3.307683000,  0.000000000],
                             [ 0.000000000,  0.000000000, 23.001852000]]);
    assert_eq!(outcar.ext_pressure, vec![-0.68, -1.59, -1.61]);
    assert_eq!(outcar.ions_per_type, vec![2, 1]);
    assert_eq!(outcar.ion_types, vec!["Se", "V"]);
    assert_eq!(outcar.ion_masses, vec![78.96, 78.96, 50.941]);
    assert_eq!(outcar.ion_iters.len(), 3);
    assert_eq!(outcar.vib, None);

    outcar.ion_iters.iter()
                    .zip(vec![27i32, 6, 4].iter())
                    .for_each(|(x, y)| assert_eq!(&x.nscf, y));

    outcar.ion_iters.iter()
                    .zip(vec![-18.95794080,
                              -18.95854979,
                              -18.95862392].iter())
                    .for_each(|(x, y)| assert_eq!(&x.toten, y));

    outcar.ion_iters.iter()
                    .zip(vec![-18.95729223,
                              -18.95789288,
                              -18.95796667].iter())
                    .for_each(|(x, y)| assert_eq!(&x.toten_z, y));

    assert_eq!(outcar.ion_iters.last().unwrap()
               .positions.last().unwrap(), &[0.00000, 0.00000, 4.13794]);
    assert_eq!(outcar.ion_iters.last().unwrap()
               .forces.last().unwrap(), &[0.000000, 0.00000 , -0.000349]);

    outcar.ion_iters.iter()
                    .zip(vec![Some(vec![0.6003306]),
                              Some(vec![0.5997977]),
                              Some(vec![0.5995733])].iter())
                    .for_each(|(x, y)| assert_eq!(&x.magmom, y));

    Ok(())
}


#[test]
fn test_ncl_outcar() -> io::Result<()> {
    let fname = get_fpath_in_current_dir!("OUTCAR_ncl");
    let outcar = Outcar::from_file(&fname)?;

    assert_eq!(outcar.lsorbit, true);
    assert_eq!(outcar.ispin, 1);
    assert_eq!(outcar.ibrion, -1);
    assert_eq!(outcar.nions, 3);
    assert_eq!(outcar.nkpts, 81);
    assert_eq!(outcar.nbands, 28);
    assert_eq!(outcar.efermi, -2.2555);
    assert_eq!(outcar.cell, [[2.864537506,-1.653841500, 0.000000000],
                             [0.000000000, 3.307683000, 0.000000000],
                             [0.000000000, 0.000000000,23.001852000]]);
    assert_eq!(outcar.ext_pressure, vec![-1.77]);
    assert_eq!(outcar.ions_per_type, vec![2, 1]);
    assert_eq!(outcar.ion_types, vec!["Se", "V"]);
    assert_eq!(outcar.ion_masses, vec![78.96, 78.96, 50.941]);
    assert_eq!(outcar.ion_iters.len(), 1);
    assert_eq!(outcar.vib, None);

    outcar.ion_iters.iter()
                    .zip(vec![31i32].iter())
                    .for_each(|(x, y)| assert_eq!(&x.nscf, y));

    outcar.ion_iters.iter()
                    .zip(vec![-19.00260977].iter())
                    .for_each(|(x, y)| assert_eq!(&x.toten, y));

    outcar.ion_iters.iter()
                    .zip(vec![-19.00194579].iter())
                    .for_each(|(x, y)| assert_eq!(&x.toten_z, y));

    assert_eq!(outcar.ion_iters.last().unwrap().positions, vec![[1.90969, -0.00000, 2.55994],
                                                                [0.95485,  1.65384, 5.71537],
                                                                [2.86454, -1.65384, 4.13794]]);

    assert_eq!(outcar.ion_iters.last().unwrap().forces, vec![[-0.000003, -0.000004, -0.012396],
                                                             [-0.000003,  0.000007,  0.013014],
                                                             [ 0.000006, -0.000003, -0.000618]]);

    outcar.ion_iters.iter()
                    .zip(vec![Some(vec![ 0.0000227, -0.0001244,  0.5998908])].iter())
                    .for_each(|(x, y)| assert_eq!(&x.magmom, y));
    Ok(())
}