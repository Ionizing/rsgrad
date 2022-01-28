use std::cmp::Ordering;

use anyhow::{anyhow, Context};
use crate::traits::Result;
use crate::types::{
    Structure,
    Mat33,
    MatX3,
};

pub struct Poscar {  // I have no plan to support vasp4 format
    pub comment: String,
    pub scale: f64,
    pub cell: Mat33<f64>,
    pub ion_types: Vec<String>,
    pub ions_per_type: Vec<i32>,
    pub pos_cart: MatX3<f64>,
    pub pos_frac: MatX3<f64>,
    pub constraints: Option<MatX3<bool>>,
}


impl Poscar {
    fn from_str(txt: &str) -> Result<Self> {
        let mut lines = txt.lines();
        let comment: String = lines.next().context("[POSCAR]: File may be blank.")?.trim().to_string();
        let scale: f64 = lines.next().context("[POSCAR]: Cannot parse scale constant.")?
            .split_whitespace()
            .next().context("[POSCAR]: Scale line may be empty.")?
            .parse::<f64>()
            .context("[POSCAR]: Scale constant cannot be converted to float number.")?;
        
        let scale = match scale.partial_cmp(&0.0) {
            Some(Ordering::Greater) => scale,
            Some(Ordering::Less) | Some(Ordering::Equal) => 
                return Err(anyhow!("[POSCAR]: Scale constant should be greater than 0.0.")),
            None => return Err(anyhow!("[POSCAR]: Scale constant cannot be NaN.")),
        };
        
        let cell: Mat33<f64> = {
            let mut v = [[0.0f64; 3]; 3];
            for i in 0..3 {
                let line = lines.next().context("[POSCAR]: Incomplete lines for cell info.")?;
                let row = line.split_whitespace().take(3).collect::<Vec<_>>();
                if row.len() < 3 {
                    return Err(anyhow!("[POSCAR]: Cell lines incomplete."));
                }
                for (j, x) in row.into_iter().enumerate() {
                    let val = x.parse::<f64>().context("[POSCAR]: Cell lines contain invalid value.")?;
                    if val.is_nan() {
                        return Err(anyhow!("[POSCAR]: Cell lines contain NaN value."));
                    }
                    v[i][j] = val;
                }
            }
            v
        };
        
        let ion_types = {
            let words = lines.next()
                .context("[POSCAR]: Element tags line not found, rsgrad has no plan to support vasp4 format.")?
                .split_whitespace()
                .take_while(|x| !x.contains("!"))
                .collect::<Vec<_>>();
            if words.len() == 0 {
                return Err(anyhow!("[POSCAR]: At lease one element is needed."));
            }
            words
        };

        let ions_per_type = {
            let numbers = lines.next()
                .context("[POSCAR]: Count of each element not found.")?
                .split_whitespace()
                .take_while(|x| !x.contains("!"))
                .map(|x| x.parse::<i32>().context("[POSCAR]: Invalid atom count of element."))
                .collect::<Result<Vec<_>>>()?;
            if numbers.len() != ion_types.len() {
                return Err(anyhow!("[POSCAR]: Inconsistent element types and atom counts."));
            }
            numbers
        };

        let line = lines.next().context("[POSCAR]: Constraints or Coordination type not found.")?;
        let has_constraints = {
            match line.trim_start().chars().next() {
                Some('s') | Some('S') => true,
                Some('d') | Some('D') | Some('c') | Some('C') | Some('k') | Some('K') => false,
                _ => return Err(anyhow!("[POSCAR]: Constraints line or Coordination type line missing."))
            }
        };
        
        let is_direct = {
            let line = if has_constraints {
                lines.next().context("[POSCAR]: Coordination type not found.")?
            } else {
                line
            };
            match line.trim_start().chars().next() {
                Some('d') | Some('D') => true,
                Some('c') | Some('C') | Some('k') | Some('K') => false,
                _ => return Err(anyhow!("[POSCAR]: Coordination type line missing."))
            }
        };

        todo!()
    }
}
