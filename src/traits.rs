use anyhow::Result;

pub trait OptProcess {
    fn process(&self) -> Result<()>;
}
