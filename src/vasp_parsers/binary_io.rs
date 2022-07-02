//! Binary read (no write stuff for now) trait, produces 1D to 3D NDArray

use std::io::{
    self,
    Result,
};

use byteorder::{LittleEndian, ReadBytesExt};
use ndarray::{
    Array1,
    Array2,
    Array3,
};
use paste::paste;

macro_rules! impl_read_1darray {
    ($t: tt) => {
        paste! {
            fn [<read_array_1d_ $t>](&mut self, len: usize) -> Result<Array1<$t>> {
                let mut ret = Array1::zeros(len);
                self.[<read_ $t _into>]::<LittleEndian>(ret.as_slice_mut().unwrap())?;
                Ok(ret)
            }
        }
    };
}

macro_rules! impl_read_2darray {
    ($t: tt) => {
        paste! {
            fn [<read_array_2d_ $t>](&mut self, nrow: usize, ncol: usize) -> Result<Array2<$t>> {
                let mut ret = Array2::zeros((nrow, ncol));
                self.[<read_ $t _into>]::<LittleEndian>(ret.as_slice_mut().unwrap())?;
                Ok(ret)
            }
        }
    };
}

macro_rules! impl_read_3darray {
    ($t: tt) => {
        paste! {
            fn [<read_array_3d_ $t>](&mut self, ni: usize, nj: usize, nk: usize) -> Result<Array3<$t>> {
                let mut ret = Array3::zeros((ni, nj, nk));
                self.[<read_ $t _into>]::<LittleEndian>(ret.as_slice_mut().unwrap())?;
                Ok(ret)
            }
        }
    };
}

pub trait ReadArray: io::Read {
    impl_read_1darray!(i32);
    impl_read_1darray!(f32);
    impl_read_1darray!(i64);
    impl_read_1darray!(f64);

    impl_read_2darray!(i32);
    impl_read_2darray!(f32);
    impl_read_2darray!(i64);
    impl_read_2darray!(f64);

    impl_read_3darray!(i32);
    impl_read_3darray!(f32);
    impl_read_3darray!(i64);
    impl_read_3darray!(f64);
}

impl<R: io::Read + ?Sized> ReadArray for R {}
