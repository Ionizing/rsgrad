pub type Result<T> = anyhow::Result<T>;

pub trait OptProcess {
    fn process(&self) -> Result<()>;
}

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
