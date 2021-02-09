use std::io::Result;
use std::path::Path;
use std::time;
use clap::{Arg, App};
use rsgrad::outcar::Outcar;
use rsgrad::outcar::PrintOptIterations;

fn main() -> Result<()> {
    let now = time::Instant::now();

    let matches = App::new("rsgrad")
        .version("0.1")
        .author("Ionizing PeterSmith_9@outlook.com")
        .about("Tracking the relaxation or MD progress of VASP calculation")
        .arg(Arg::with_name("input")
             .default_value("OUTCAR"))
        .get_matches();

    let f = Outcar::from_file(Path::new(matches.value_of("input").unwrap()))?;
    println!("{}", PrintOptIterations::from(f.ion_iters));

    eprintln!("# Time used: {:?}", now.elapsed());
    Ok(())
}
