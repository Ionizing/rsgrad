use std::time;

use env_logger::init_from_env;
use log::info;
use rsgrad::{
    Result,
    cli,
};



fn main() -> Result<()> {
    let now = time::Instant::now();

    init_from_env(
        env_logger::Env::new().filter_or("RSGRAD_LOG", "info"));

    cli::run()?;

    info!("Time used: {:?}", now.elapsed());
    Ok(())
}
