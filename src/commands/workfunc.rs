use std::path::PathBuf;
use clap::Args;
use rayon;
use anyhow::{
    Context,
    anyhow,
};
use log::info;
use ndarray;
use plotly;

use crate::{
    types::{
        Result,
        Axis,
    },
    OptProcess,
    ChargeDensity,
    ChargeType,
    Outcar,
    commands::common::write_array_to_txt,
};

#[derive(Debug, Args)]
/// Calculate work-function from LOCPOT file, OUTCAR is also needed to get the Fermi level.
///
/// The work function is calculated by plannar integration of the data cube in LOCPOT. The
/// selected axis should be perpendicular to the other two axises.
pub struct Workfunc {
    #[arg(default_value="./LOCPOT")]
    /// LOCPOT file path. Turn on 'LVHAR' in INCAR to get the electro-static potential saved it.
    locpot: PathBuf,

    #[arg(long, short = 'o', default_value="./workfunction.html")]
    /// Write the plot to html and view it in the web browser.
    htmlout: PathBuf,

    #[arg(long, default_value="locpot.txt")]
    /// Write the raw plot data as txt file in order to replot it with more advanced tools.
    txtout: PathBuf,

    #[arg(long, default_value="./OUTCAR")]
    /// OUTCAR file path. This file is needed to get the E-fermi level and lattice properties.
    outcar: PathBuf,

    #[arg(long, default_value="z", value_enum, ignore_case = true)]
    /// Integration direction. e.g. if 'z' is provided, the XoY plane is integrated.
    axis: Axis,

    #[arg(long)]
    /// Open default browser to see the plot immediately.
    show: bool,

    #[arg(long)]
    /// Render the plot and print the rendered code to stdout.
    to_inline_html: bool,
}


impl OptProcess for Workfunc {
    fn process(&self) -> Result<()> {
        let mut locpot: Result<ChargeDensity> = Err(anyhow!(""));
        let mut outcar: Result<Outcar> = Err(anyhow!(""));

        rayon::scope(|s| {
            s.spawn(|_| {
                info!("Reading electro-static potential data from {:?}", &self.locpot);
                locpot = ChargeDensity::from_file(&self.locpot, ChargeType::Locpot);
            });
            s.spawn(|_| {
                info!("Reading {:?}", &self.outcar);
                outcar = Outcar::from_file(&self.outcar);
            });
        });

        let locpot = locpot.context(format!("Parse file {:?} failed.", self.locpot))?;
        let outcar = outcar.context(format!("Parse file {:?} failed.", self.outcar))?;

        let efermi = outcar.efermi;
        let ngrid = locpot.ngrid;
        let cell = locpot.pos.cell;
        let iaxis = match self.axis {
            Axis::X => 0usize,
            Axis::Y => 1usize,
            Axis::Z => 2usize,
        };
        let axislen = {
            let row = cell[iaxis];
            (row[0] * row[0] + row[1] * row[1] + row[2] * row[2]).sqrt()
        };

        let workfunc = match self.axis {
            Axis::X => {
                locpot.chg[0]
                    .mean_axis(ndarray::Axis(2)).unwrap()
                    .mean_axis(ndarray::Axis(1)).unwrap()
            },
            Axis::Y => {
                locpot.chg[0]
                    .mean_axis(ndarray::Axis(2)).unwrap()
                    .mean_axis(ndarray::Axis(0)).unwrap()
            },
            Axis::Z => {
                locpot.chg[0]
                    .mean_axis(ndarray::Axis(1)).unwrap()
                    .mean_axis(ndarray::Axis(0)).unwrap()
            },
        } - efermi;

        let distance = ndarray::Array::linspace(0.0, axislen, ngrid[iaxis]);

        info!("Writing raw plot data to {:?}", self.txtout);
        write_array_to_txt(&self.txtout, vec![&distance, &workfunc], "Distance(A)  E-Ef(eV)")?;

        let trace = plotly::Scatter::from_array(distance, workfunc)
            .mode(plotly::common::Mode::Lines);

        let mut plot = plotly::Plot::new();
        plot.add_trace(trace);
        plot.use_local_plotly();

        let layout = plotly::Layout::new()
            .title(plotly::common::Title::new(&format!("Work function along {} axis", self.axis)))
            .y_axis(plotly::layout::Axis::new()
                    .title(plotly::common::Title::new("E-Ef (eV)"))
                    .zero_line(true))
            .x_axis(plotly::layout::Axis::new()
                    .title(plotly::common::Title::new("Distance (A)"))
                    .zero_line(true));
        plot.set_layout(layout);

        info!("Writing to {:?}", self.htmlout);
        plot.to_html(&self.htmlout);

        if self.show {
            plot.show();
        }

        if self.to_inline_html {
            info!("Printing inline html to stdout ...");
            println!("{}", plot.to_inline_html(None));
        }

        Ok(())
    }
}
