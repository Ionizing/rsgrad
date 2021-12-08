use crate::types::{
    Matrix,
    Vector,
};

#[derive(Clone)]
pub struct Kpoints {
    pub nkpoints:   u32,
    pub kpointlist: Matrix<f64>,
    pub weights:    Vector<f64>,
}

