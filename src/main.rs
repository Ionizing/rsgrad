use std::io::Result;
use std::path::Path;
use std::time;
use clap::{Arg, App, AppSettings};
use rsgrad::outcar::Outcar;
use rsgrad::format::IonicIterationsFormat;


fn main() -> Result<()> {
    let now = time::Instant::now();

    let matches = App::new("rsgrad")
        .setting(AppSettings::ColoredHelp)
        .version("0.2")
        .author("Ionizing PeterSmith_9@outlook.com")
        .about("Tracking the relaxation or MD progress of VASP calculation")
        .arg(Arg::with_name("input")
             .default_value("OUTCAR"))
        .arg(Arg::from_usage("--print-energy [t|f] 'Print TOTEN in eV'")
             .default_value("f"))
        .arg(Arg::from_usage("--print-energyz [t|f] 'Print TOTEN without entropy in eV'")
             .default_value("t"))
        .arg(Arg::from_usage("--print-lgde [t|f] 'Print Log10(deltaE), where deltaE is the
    difference between two consecutive energies
    without entropy'")
             .default_value("t"))
        .arg(Arg::from_usage("--print-favg [t|f] 'Print averaged total force in eV/A'")
             .default_value("t"))
        .arg(Arg::from_usage("--print-fmax [t|f] 'Print maximum total force in eV/A'")
             .default_value("t"))
        .arg(Arg::from_usage("--print-fmax-index [t|f] 'Print the index of ion with maximum force load'")
             .default_value("f"))
        .arg(Arg::from_usage("--print-fmax-axis [t|f] 'Print the axis where the strongest force weight lies on'")
             .default_value("f"))
        .arg(Arg::from_usage("--print-nscf [t|f] 'Print the number SCF iteration in the ionic iteration'")
             .default_value("t"))
        .arg(Arg::from_usage("--print-time [t|f] 'Print the time usage of each ionic iteration in unit of minute'")
             .default_value("t"))
        .arg(Arg::from_usage("--print-magmom [t|f] 'Print the magnetic moment in current system'")
             .default_value("t"))
        .arg(Arg::from_usage("--print-volume [t|f] 'Print the system volumes of each ionic iteration'")
             .default_value("f"))
        .get_matches();

    let convert_helper = |c: &str| -> bool {
        match c {
            "t" | "T" => true,
            "f" | "F" => false,
            _ => unreachable!("Only [tTfF] is allowed in --print-xxx values")
        }
    };

    let f = Outcar::from_file(Path::new(matches.value_of("input").unwrap()))?;

    let iif = IonicIterationsFormat::from(f.ion_iters)
        .print_energy     (convert_helper(matches.value_of("print-energy").unwrap()))
        .print_energyz    (convert_helper(matches.value_of("print-energyz").unwrap()))
        .print_log10de    (convert_helper(matches.value_of("print-lgde").unwrap()))
        .print_favg       (convert_helper(matches.value_of("print-favg").unwrap()))
        .print_fmax       (convert_helper(matches.value_of("print-fmax").unwrap()))
        .print_fmax_axis  (convert_helper(matches.value_of("print-fmax-axis").unwrap()))
        .print_fmax_index (convert_helper(matches.value_of("print-fmax-index").unwrap()))
        .print_nscf       (convert_helper(matches.value_of("print-nscf").unwrap()))
        .print_time_usage (convert_helper(matches.value_of("print-time").unwrap()))
        .print_magmom     (convert_helper(matches.value_of("print-magmom").unwrap()))
        .print_volume     (convert_helper(matches.value_of("print-volume").unwrap()));

    print!("{}", iif);

    eprintln!("# Time used: {:?}", now.elapsed());
    Ok(())
}
