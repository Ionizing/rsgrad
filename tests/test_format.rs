use std::path::PathBuf;
use std::io;
use std::fs;
use rsgrad::{
    outcar::Outcar,
    format::Trajectory,
    format::Vibrations,
};

use vasp_poscar::Poscar;
use tempdir::TempDir;

macro_rules! get_fpath_in_current_dir {
    ($fname:expr) => {{
        let mut path = PathBuf::from(file!());
        path.pop();
        path.push($fname);
        path
    }}
}

#[test]
fn test_save_as_xdatcar() -> io::Result<()> {
    let fname = get_fpath_in_current_dir!("OUTCAR_another_rlx");
    let outcar = Outcar::from_file(&fname)?;
    let traj = Trajectory::from(outcar);

    let tmpdir = TempDir::new("rsgrad_test")?;
    let path = tmpdir.path();
    traj.save_as_xdatcar(&path)?;

    let xdatcar_ref = fs::read_to_string(
        get_fpath_in_current_dir!("XDATCAR_another_rlx"))?;
    let xdatcar_content = fs::read_to_string(path.join("XDATCAR"))?;
    assert_eq!(xdatcar_content, xdatcar_ref);
    Ok(())
}

#[test]
fn test_save_as_seperated_poscars() -> io::Result<()> {
    let fname = get_fpath_in_current_dir!("OUTCAR_another_rlx");
    let outcar = Outcar::from_file(&fname)?;
    let traj = Trajectory::from(outcar);

    let tmpdir = TempDir::new("rsgrad_test")?;
    traj.save_as_seperated_poscars(tmpdir.path())?;

    // Validation
    let entries = fs::read_dir(tmpdir)?
        .map(|res| res.map(|e| e.path()))
        .collect::<Result<Vec<_>, io::Error>>()?;

    assert!(entries.iter().all(|f| Poscar::from_path(f).is_ok()));
    Ok(())
}

#[test]
fn test_save_as_single_xsf() -> io::Result<()> {
    let fname = get_fpath_in_current_dir!("OUTCAR_vibrations");
    let outcar = Outcar::from_file(&fname)?;
    let vibs = Vibrations::from(outcar);

    let tmpdir = TempDir::new_in(".", "rsgrad_test").unwrap();
    vibs.save_as_xsf(1, &tmpdir.path())?;

    for i in 1..=vibs.modes.len() {
        vibs.save_as_xsf(i, &tmpdir.path())?;
    }

    Ok(())
}
