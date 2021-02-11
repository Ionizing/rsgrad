use std::io::Result;
use std::path::Path;
use std::time;
use clap::{Arg, App};
use rsgrad::outcar::Outcar;
// use rsgrad::outcar::PrintOptIterations;
use rsgrad::format::IonicIterationsFormat;


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

    let iif = IonicIterationsFormat::from(f.ion_iters)
        .print_energy(true)
        .print_energyz(true)
        .print_log10de(true)
        .print_favg(true)
        .print_fmax(true)
        .print_fmax_axis(true)
        .print_fmax_index(true)
        .print_nscf(true)
        .print_time_usage(true)
        .print_magmom(true)
        .print_volume(true);

    print!("{}", iif);

    eprintln!("# Time used: {:?}", now.elapsed());
    Ok(())
}
