use std::io::Result;
use std::path::Path;
use std::time;
use clap::{Arg, App};
use rsgrad::outcar::Outcar;

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
    for (i, ionit) in f.ion_iters.iter().enumerate() {
        println!("{:3} {}", i+1, ionit);
    }

    eprintln!("# Time used: {:?}", now.elapsed());
    Ok(())
}
