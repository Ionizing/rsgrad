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
