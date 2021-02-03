use std::io::Result;
use std::path::Path;
use rsgrad::outcar::Outcar;

fn main() -> Result<()> {
    let f = Outcar::from_file(Path::new("OUTCAR"))?;

    Ok(())
}
