use ndarray::{
    Array1,
    Array2,
    Array3,
};

pub type Vector<T> = Array1<T>;  // Define this type to use broadcast operations.
pub type Matrix<T> = Array2<T>;
pub type Cube<T>   = Array3<T>;
pub type MatX3<T> = Vec<[T;3]>;  // Nx3 matrix
pub type Mat33<T> = [[T;3];3];   // 3x3 matrix


#[derive(Clone)]
pub struct Structure {
    pub cell          : Mat33<f64>,
    pub ion_types     : Vec<String>,
    pub ions_per_type : Vec<i32>,
    pub car_pos       : MatX3<f64>,
    pub frac_pos      : MatX3<f64>,
    pub constr        : Option<MatX3<bool>>,
}
