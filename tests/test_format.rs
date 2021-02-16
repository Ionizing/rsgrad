use std::path::PathBuf;
use std::io;
use std::fs;
use rsgrad::{
    outcar::Outcar,
    format::Trajectory,
};

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
    let path = tmpdir.path().join("XDATCAR");
    traj.save_as_xdatcar(&path)?;

    let xdatcar_ref = fs::read_to_string(
        get_fpath_in_current_dir!("XDATCAR_another_rlx"))?;
    let xdatcar_content = fs::read_to_string(path)?;
    assert_eq!(xdatcar_content, xdatcar_ref);
    Ok(())
}
