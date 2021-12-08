use ndarray::Array5;
use crate::types::{
    //Vector,
    //Matrix,
    Cube,
};

#[derive(Clone)]
pub struct ProjectedDOS {
    pub nions:       u32,
    pub nspin:       u32,
    pub nkpoints:    u32,
    pub nbands:      u32,
    pub lsorbit:     bool,
    pub nlm:         Vec<String>,
    pub eigvals:     Cube<f64>,             // [ispin, ikpoint, iband]
    pub occupations: Cube<f64>,
    pub projected:   Array5<f64>,     // [[ispin, ikpoint, iband], [iion, iorbit]]
}
