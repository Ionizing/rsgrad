use std::path::PathBuf;
use structopt::{
    StructOpt,
    clap::AppSettings,
};
use log::{
    warn, 
    info,
    error,
};
use anyhow::bail;
use crate::{
    Poscar,
    Result,
    OptProcess,
    index_transform,
};


#[derive(Debug, StructOpt)]
#[structopt(setting = AppSettings::ColoredHelp,
            setting = AppSettings::ColorAuto,
            setting = AppSettings::AllowNegativeNumbers)]
/// Operation(s) about POSCAR, including split it into two POSCARs.
pub struct Pos {
    #[structopt(default_value = "./POSCAR")]
    /// Specify the input POSCAR file
    poscar: PathBuf,

    #[structopt(short = "i", long)]
    /// Selects the indices to operate.
    ///
    /// Step indices start from '1', if '0' is given, all the structures will be selected.
    /// Step indices can be negative, where negative index means counting reversely.
    /// E.g. "-i -2 -1 1 2 3" means selecting the last two and the first three atom.
    select_indices: Option<Vec<i32>>,

    #[structopt(short = "a", long, default_value = "POSCAR_A")]
    /// Splitted POSCAR path with selected atoms
    a_name: PathBuf,

    #[structopt(short = "b", long, default_value = "POSCAR_B")]
    /// Splitted POSCAR path with complement of `a_name`
    b_name: PathBuf,

    #[structopt(long)]
    /// The symbols of each atom will not be written as comment in POSCAR
    no_add_symbols_tags: bool,

    #[structopt(long)]
    /// Atom constraints will be dropped when writting POSCAR
    no_preserve_constraints: bool,

    #[structopt(short = "c", long)]
    /// Cartesian coordinates is used in writting POSCAR
    cartesian: bool,

    #[structopt(short = "s", long)]
    /// Split POSCAR according to selected_indices
    split: bool,

    #[structopt(long)]
    /// Convert POSCAR to cartesian coordinates or fractional coordinates
    convert: bool,

    #[structopt(long, default_value = "POSCAR_new")]
    /// The target path of converted POSCAR
    converted: PathBuf,
}


fn poscar_split(poscar: &Poscar, inds: &[usize]) -> (Poscar, Poscar) {
    // Assume the inds is neither empty nor full.
    
    let inds_a = inds.iter()
        .map(|x| *x as usize)
        .collect::<Vec<_>>();
    let inds_b = (0 .. poscar.pos_cart.len())
        .into_iter()
        .filter(|i| !inds.contains(i))
        .collect::<Vec<_>>();

    let (comment_a, comment_b) = ("Generated by rsgrad, POSCAR with selected atoms",
                                  "Generated by rsgrad, POSCAR complement");
    let scale = poscar.scale;
    let cell = poscar.cell;

    let ion_types_total = &poscar.ion_types;
    let ions_per_type_total = &poscar.ions_per_type;

    let cumsum = poscar.ions_per_type.iter()
        .scan(0, |acc, x| {
            *acc += *x as usize;
            Some(*acc)
        })
        .collect::<Vec<usize>>();

    let (ion_types_a, ions_per_type_a, ion_types_b, ions_per_type_b) = {
        let mut ions_per_type = vec![0usize; ions_per_type_total.len()];
        for _i in inds_a.iter() {
            let i = cumsum.binary_search(&(_i + 1)).unwrap_or_else(|x| x);
            ions_per_type[i] += 1;
        }
        let ion_types_a = ions_per_type.iter().zip(ion_types_total.iter())
            .filter(|(n, _)| n != &&0)
            .map(|(_, s)| s.clone())
            .collect::<Vec<_>>();
        let ions_per_type_a = ions_per_type.iter()
            .filter(|x| x > &&0)
            .map(|x| *x as i32)
            .collect::<Vec<i32>>();

        let ions_per_type_b = ions_per_type_total
            .iter().zip(ions_per_type.iter())
            .map(|(t, a)| *t - *a as i32)
            .collect::<Vec<i32>>();

        let ion_types_b = ions_per_type_b.iter().zip(ion_types_total.iter())
            .filter(|(n, _)| n != &&0)
            .map(|(_, s)| s.clone())
            .collect::<Vec<_>>();

        let ions_per_type_b = ions_per_type_b.into_iter()
            .filter(|x| x > &0)
            .collect::<Vec<i32>>();

        (ion_types_a, ions_per_type_a, ion_types_b, ions_per_type_b)
    };
        
    let pos_cart_a = inds_a.iter()
        .map(|i| poscar.pos_cart[*i])
        .collect::<Vec<_>>();

    let pos_frac_a = inds_a.iter()
        .map(|i| poscar.pos_frac[*i])
        .collect::<Vec<_>>();
    
    let pos_cart_b = inds_b.iter()
        .map(|i| poscar.pos_cart[*i])
        .collect::<Vec<_>>();

    let pos_frac_b = inds_b.iter()
        .map(|i| poscar.pos_frac[*i])
        .collect::<Vec<_>>();

    let (constraints_a, constraints_b) = if let Some(constraints) = poscar.constraints.as_ref() {
        (
            Some(inds_a.iter()
                    .map(|i| constraints[*i])
                    .collect::<Vec<_>>())
            ,
            Some(inds_b.iter()
                    .map(|i| constraints[*i])
                    .collect::<Vec<_>>())
        )
    } else {
        (None, None)
    };

    let poscar_a = Poscar {
        comment: comment_a.into(),
        scale,
        cell,
        ion_types: ion_types_a,
        ions_per_type: ions_per_type_a,
        pos_cart: pos_cart_a,
        pos_frac: pos_frac_a,
        constraints: constraints_a,
    };

    let poscar_b = Poscar {
        comment: comment_b.into(),
        scale,
        cell,
        ion_types: ion_types_b,
        ions_per_type: ions_per_type_b,
        pos_cart: pos_cart_b,
        pos_frac: pos_frac_b,
        constraints: constraints_b,
    };

    (poscar_a, poscar_b)
}


impl OptProcess for Pos {
    fn process(&self) -> Result<()> {
        info!("Reading POSCAR file {:?} ...", &self.poscar);
        let pos = Poscar::from_file(&self.poscar)?;


        if self.convert {
            info!("Converting it to {:?}", &self.converted);

            pos.to_formatter()
                .preserve_constraints(!self.no_preserve_constraints)
                .fraction_coordinates(!self.cartesian)
                .add_symbol_tags(!self.no_add_symbols_tags)
                .to_file(&self.converted)?;
            
            info!("Done");
        }


        if self.split {
            info!("Splitting it to {:?} and {:?} ...", &self.a_name, &self.b_name);
            if let Some(select_indices) = self.select_indices.as_ref() {
                let inds = index_transform(select_indices.clone(), pos.get_natoms() as usize)
                    .into_iter()
                    .map(|i| i - 1)
                    .collect::<Vec<_>>();
                if inds == (0 .. pos.get_natoms() as usize).into_iter().collect::<Vec<_>>() {
                    error!("You selected all the atoms for {:?}, leaving {:?} nothing,\
please check you input (`0` means selecting all the atoms)",
                        &self.a_name, &self.b_name);
                    bail!("Aborted due to invalid selected indices.");
                }

                let (poscar_a, poscar_b) = poscar_split(&pos, &inds);

                info!("{:?} contains", &self.a_name);
                for (s, c) in poscar_a.ion_types.iter().zip(poscar_a.ions_per_type.iter()) {
                    info!("  {:>5}  {:>4}", s, c);
                }

                info!("{:?} contains", &self.b_name);
                for (s, c) in poscar_b.ion_types.iter().zip(poscar_b.ions_per_type.iter()) {
                    info!("  {:>5}  {:>4}", s, c);
                }

                poscar_a.to_formatter()
                    .preserve_constraints(!self.no_preserve_constraints)
                    .fraction_coordinates(!self.cartesian)
                    .add_symbol_tags(!self.no_add_symbols_tags)
                    .to_file(&self.a_name)?;
                info!("{:?} written", &self.a_name);

                poscar_b.to_formatter()
                    .preserve_constraints(!self.no_preserve_constraints)
                    .fraction_coordinates(!self.cartesian)
                    .add_symbol_tags(!self.no_add_symbols_tags)
                    .to_file(&self.b_name)?;
                info!("{:?} written", &self.b_name);
            } else {
                warn!("No atom indices selected, nothing done!");
                bail!("Aborted due to invalid selected indices.");
            }
        }

        Ok(())
    }
}
