use std::path::PathBuf;

use rsgrad::{
    Result,
    Poscar,
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
fn test_read_poscar() -> Result<()> {
    // fname volume [lattice length] [angles]
    let fnames = [
        ("POSCAR.Al12O18",              260.8221301178378, 
            [  4.79489982,   4.79489982,  13.0995], [90., 90., 119.99999751]),
        ("POSCAR.CdS_HSE",              102.15864298587704,
            [  4.17045897,   4.17045907,   6.7822969], [90., 90., 120.00000695]),
        ("POSCAR.Fe3O4",                144.57668906044427,
            [  2.95442   ,   5.11432539,   9.994342], [90., 90., 106.78838383]),
        ("POSCAR.Li2O",                 790.4268189360001,
            [ 9.246,  9.246,  9.246], [ 90.   , 90.   , 90.]),
        ("POSCAR.LiFePO4",              300.12708022907657,
            [10.41015404,  6.06327401,  4.75489403], [89.99235315, 90.0097851 , 89.99856666]),
        ("POSCAR.lobster.nonspin_DOS",   39.892877896074175,
            [ 3.83533651,  3.83533618,  3.835336], [60.00000156, 60.00000442, 60.00000506]),
        ("POSCAR.lobster.spin_DOS",      39.892877896074175,
            [ 3.83533651,  3.83533618,  3.835336], [60.00000156, 60.00000442, 60.00000506]),
        ("POSCAR.O2",                   117.70356181904295,
            [  4.17727252,   5.67453122,   5.67453122], [76.94621456, 110.12319067, 110.12319067]),
        ("POSCAR.tricky_symmetry",      301.92747145565284,
            [ 2.97098261,  2.97098176, 34.27053249], [ 92.4856315 , 92.48554312, 89.9971654]),
    ];

    let tmpdir = TempDir::new("poscar_test")?.into_path();

    for (f, v, l, a) in fnames {
        println!("testing {}", f);
        let pos = Poscar::from_file(&get_fpath_in_current_dir!(f))?;
        assert_eq!(pos.get_volume(), v);
        let (len, ang) = pos.get_cell_params();
        const THRESHOLD: f64 = 1E-6;

        assert!((len[0] - l[0]).abs() < THRESHOLD);
        assert!((len[1] - l[1]).abs() < THRESHOLD);
        assert!((len[2] - l[2]).abs() < THRESHOLD);

        assert!((ang[0] - a[0]).abs() < THRESHOLD);
        assert!((ang[1] - a[1]).abs() < THRESHOLD);
        assert!((ang[2] - a[2]).abs() < THRESHOLD);

        if vec!["POSCAR.O2"].contains(&f) {
            assert!(pos.constraints.is_some());
            let constraints = pos.constraints.as_ref().unwrap();
            assert_eq!(constraints[0], [true, true, false]);
            assert_eq!(constraints[4], [true, false, true]);
            assert_eq!(constraints[7], [false, true, true]);
        } else {
            assert!(pos.constraints.is_none());
        }

        let mut path = tmpdir.clone();
        path.push(f);
        println!("{:?}", &path);
        pos.to_formatter()
            .preserve_constraints(true)
            .fraction_coordinates(true)
            .add_symbol_tags(false)
            .to_file(&path)?;
        Poscar::from_file(&path)?;
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
