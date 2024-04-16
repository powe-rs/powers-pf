use crate::pfopt::MPOpt;
use powers::SBus;

use anyhow::{format_err, Result};
use num_complex::Complex64;
use sparsetools::csr::CSR;
use spsolve::Solver;

pub trait ProgressMonitor {
    fn update(&self, i: usize, norm_f: f64);
}

pub struct PrintProgress {}

impl ProgressMonitor for PrintProgress {
    fn update(&self, i: usize, norm_f: f64) {
        if i == 0 {
            println!(" it    max P & Q mismatch (p.u.)");
            println!("----  ---------------------------");
            // println!("%3d        %10.3e", i, normF);
            println!("{}        {}", i, norm_f);
        } else {
            // println!("%3d        %10.3e", i, normF);
            println!("{}        {}", i, norm_f);
        }
    }
}

pub(crate) fn gausspf(
    _y_bus: &CSR<usize, Complex64>,
    _s_bus: &dyn SBus,
    _v0: &[Complex64],
    _ref: &[usize],
    _pv: &[usize],
    _pq: &[usize],
    _solver: &dyn Solver<usize, f64>,
    _mpopt: &MPOpt,
    _progress: Option<&dyn ProgressMonitor>,
) -> Result<(Vec<Complex64>, bool, usize)> {
    Err(format_err!("not implemented"))
}
