use crate::pfopt::{Alg, MPOpt};
use powers::{make_ybus, SBus};

use anyhow::{format_err, Result};
use caseformat::{Branch, Bus};
use num_complex::Complex64;
use sparsetools::coo::{CCoo, Coo};
use sparsetools::csr::CSR;
use spsolve::Solver;

/// Solves the power flow using a fast decoupled method.
///
/// Solves for bus voltages given the full system admittance matrix (for
/// all buses), the complex bus power injection vector (for all buses),
/// the initial vector of complex bus voltages, the FDPF matrices B prime
/// and B double prime, and column vectors with the lists of bus indices
/// for the swing bus, PV buses, and PQ buses, respectively. The bus voltage
/// vector contains the set point for generator (including ref bus)
/// buses, and the reference angle of the swing bus, as well as an initial
/// guess for remaining magnitudes and angles. MPOPT is a MATPOWER options
/// vector which can be used to set the termination tolerance, maximum
/// number of iterations, and output options (see MPOPTION for details).
/// Uses default options if this parameter is not given. Returns the
/// final complex voltages, a flag which indicates whether it converged
/// or not, and the number of iterations performed.
pub(crate) fn fdpf(
    _y_bus: &CSR<usize, Complex64>,
    _s_bus: &dyn SBus,
    _v0: &[Complex64],
    _b_p: &CSR<usize, f64>,
    _b_pp: &CSR<usize, f64>,
    _ref: &[usize],
    _pv: &[usize],
    _pq: &[usize],
    _solver: &dyn Solver<usize, f64>,
    _mpopt: &MPOpt,
    _progress: Option<&dyn ProgressMonitor>,
) -> Result<(Vec<Complex64>, bool, usize)> {
    Err(format_err!("not implemented"))
}

/// Builds the two matrices B prime and B double prime used in the fast
/// decoupled power flow.
pub(crate) fn make_b(
    // mpc: &MPC,
    base_mva: f64,
    bus: &[Bus],
    branch: &[Branch],
    alg: Alg,
    double_prime: bool,
) -> (Coo<usize, f64>, Option<Coo<usize, f64>>) {
    // Form Bp (B prime).
    let mut bus = bus.to_vec(); // modify a copy of bus
    for b in bus.iter_mut() {
        b.bs = 0.0; // zero out shunts at buses
    }

    let b_p = {
        let mut branch = branch.to_vec(); // modify a copy of branch
        for br in branch.iter_mut() {
            br.br_b = 0.0; // zero out line charging shunts
            br.tap = 1.0; // cancel out taps
            if alg == Alg::FDXB {
                br.br_r = 0.0; // zero out line resistance
            }
        }
        let (y_p, _, _) = make_ybus(base_mva, &bus, &branch, false);
        // y_p.map(|y| -y.im)
        -y_p.imag()
    };

    let b_pp = if double_prime {
        // Form Bpp (B double prime).
        let mut branch = branch.to_vec(); // modify a copy of branch
        for br in branch.iter_mut() {
            br.shift = 0.0; // zero out phase shifters
            if alg == Alg::FDBX {
                br.br_r = 0.0; // zero out line resistance
            }
        }
        let (y_pp, _, _) = make_ybus(base_mva, &bus, &branch, false);
        // Some(y_pp.map(|y| -y.im))
        Some(-y_pp.imag())
    } else {
        None
    };

    (b_p, b_pp)
}

pub trait ProgressMonitor {
    fn update(&self, i: usize, norm_p: f64, norm_q: f64, p_update: bool);
}

pub struct PrintProgress {}

impl ProgressMonitor for PrintProgress {
    fn update(&self, i: usize, norm_p: f64, norm_q: f64, p_update: bool) {
        if i == 0 {
            println!("iteration     max mismatch (p.u.)  ");
            println!("type   #        P            Q     ");
            println!("---- ----  -----------  -----------");
            // println!("  -  %3d   %10.3e   %10.3e", i, norm_p, norm_q);
            println!("  -  {}   {}   {}", i, norm_p, norm_q);
        } else {
            if p_update {
                // println!("  P  %3d   %10.3e   %10.3e", i, norm_p, norm_q);
                println!("  P  {}   {}   {}", i, norm_p, norm_q);
            } else {
                // println!("  Q  %3d   %10.3e   %10.3e", i, norm_p, norm_q);
                println!("  Q  {}   {}   {}", i, norm_p, norm_q);
            }
        }
    }
}
