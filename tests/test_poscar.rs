use std::path::PathBuf;

use rsgrad::traits::Result;
use rsgrad::vasp_parsers::poscar::Poscar;


macro_rules! get_fpath_in_current_dir {
    ($fname:expr) => {{
        let mut path = PathBuf::from(file!());
        path.pop();
        path.push($fname);
        path
    }}
}


#[test]
fn test_read_poscar() -> Result<()> {
    let fnames = [
        "POSCAR.Al12O18",
        "POSCAR.CdS_HSE",
        "POSCAR.Fe3O4",
        "POSCAR.Li2O",
        "POSCAR.LiFePO4",
        "POSCAR.lobster.nonspin_DOS",
        "POSCAR.lobster.spin_DOS",
        "POSCAR.O2",
        "POSCAR.tricky_symmetry",
    ];
    for f in fnames {
        println!("testing {}", f);
        Poscar::from_file(&get_fpath_in_current_dir!(f))?;
    }

    Ok(())
}


#[test]
#[should_panic]
fn test_read_failed() {
    let fnames = [
        "POSCAR",
        "POSCAR.symbols_natoms_multilines"
    ];
    for f in fnames {
        println!("testing {}", f);
        Poscar::from_file(&get_fpath_in_current_dir!(f)).unwrap();
    }
}
