use ndarray::{
    Array1,
    Array2,
    Array3,
};

pub type Result<T> = anyhow::Result<T>;

pub trait OptProcess {
    fn process(&self) -> Result<()>;
}

/// Index array containing negative indices => Index array full of positive indices.
/// `-1` means the last index,
/// If `v` contains `0`, selecting the total indices,  `1..=len` is returned.
pub fn index_transform(v: Vec<i32>, len: usize) -> Vec<usize> {
    if v.contains(&0) {
        (1 ..= len).collect()
    } else {
        v.into_iter()
         .map(|i| {
            if i < 0 {
                i.rem_euclid(len as i32) as usize + 1
            } else {
                i as usize
            }
         })
        .collect()
    }
}

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
