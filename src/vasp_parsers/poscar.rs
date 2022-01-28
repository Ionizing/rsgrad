use crate::types::{
    Structure,
    Mat33,
    MatX3,
};

pub struct Poscar {
    comment: String,
    scale: f64,
    cell: Mat33<f64>,
    ion_types: Vec<String>,
    ions_per_type: Vec<i32>,
    pos_cart: MatX3<f64>,
    pos_frac: MatX3<f64>,
    constraints: Option<MatX3<bool>>,
}
