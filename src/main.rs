use std::io::Result;
use std::path::Path;
use rsgrad::outcar::Outcar;

fn main() -> Result<()> {
    let f = Outcar::from_file(Path::new("OUTCAR_multiple_ionic_steps"))?;
    for ionit in f.ion_iters {
        println!("{}", ionit);
    }
    Ok(())
}
