use std::iter::zip;

use anyhow::Result;
use caseformat::{Branch, Bus};
use full::iter::neg;
use full::slice::{abs, angle, norm_inf, polar};
use itertools::izip;
use num_complex::Complex64;
use powers::debug::format_rect_vec;
use powers::{make_ybus, SBus};
use sparsetools::csr::CCSR;
use sparsetools::csr::CSR;
use spsolve::FactorSolver;

use crate::pfopt::{Alg, MPOpt};

/// Solves the power flow using a fast decoupled method.
///
/// Solves for bus voltages given the full system admittance matrix (for
/// all buses), the complex bus power injection vector (for all buses),
/// the initial vector of complex bus voltages, the FDPF matrices `B` prime
/// and `B` double prime, and column vectors with the lists of bus indices
/// for the swing bus, PV buses, and PQ buses, respectively. The bus voltage
/// vector contains the set point for generator (including ref bus)
/// buses, and the reference angle of the swing bus, as well as an initial
/// guess for remaining magnitudes and angles. `mpopt` is a MATPOWER options
/// vector which can be used to set the termination tolerance, maximum
/// number of iterations, and output options.
/// Uses default options if this parameter is not given. Returns the
/// final complex voltages, a flag which indicates whether it converged
/// or not, and the number of iterations performed.
#[allow(non_snake_case)]
pub(crate) fn fdpf<F>(
    Ybus: &CSR<usize, Complex64>,
    Sbus: &dyn SBus,
    V0: &[Complex64],
    Bp: &CSR<usize, f64>,
    Bpp: &CSR<usize, f64>,
    _ref: &[usize],
    pv: &[usize],
    pq: &[usize],
    solver: &dyn FactorSolver<usize, f64, F>,
    mpopt: &MPOpt,
    progress: Option<&dyn ProgressMonitor>,
) -> Result<(Vec<Complex64>, bool, usize)> {
    let pv_pq = [pv, pq].concat();

    // options
    let tol = mpopt.pf.tolerance;
    let max_it = mpopt.pf.max_it_fd;

    // initialize
    let mut converged = false;
    let mut i = 0;
    let mut V = V0.to_vec();
    let mut Va = angle(&V);
    let mut Vm = abs(&V);

    // evaluate initial mismatch
    let (mut P, mut Q) = {
        let Ibus: Vec<Complex64> = Ybus * &V;
        let Sbus = Sbus.s_bus(&Vm);
        log::trace!("Sbus0: {}", format_rect_vec(&Sbus));
        // mis = (V .* conj(Ybus * V) - Sbus(Vm)) ./ Vm;
        let mis: Vec<Complex64> = izip!(&V, &Ibus, &Sbus, &Vm)
            .map(|(V, Ibus, Sbus, Vm)| V * Ibus.conj() - Sbus / Vm)
            .collect();
        (
            pv_pq.iter().map(|&i| mis[i].re).collect::<Vec<_>>(),
            pq.iter().map(|&i| mis[i].im).collect::<Vec<_>>(),
        )
    };

    // check tolerance
    let normP = norm_inf(&P);
    let normQ = norm_inf(&Q);
    if let Some(pm) = progress {
        pm.update(i, i, normP, normQ);
    }
    if normP < tol && normQ < tol {
        converged = true;
        log::info!("Converged!");
    }
    log::debug!("normP = {}, normQ = {}", normP, normQ);

    // reduce B matrices
    let Bp = Bp.select(Some(&pv_pq), Some(&pv_pq))?;
    let Bpp = Bpp.select(Some(&pq), Some(&pq))?;

    let LUp = solver.factor(Bp.cols(), Bp.colidx(), Bp.rowptr(), Bp.values())?;
    let LUpp = solver.factor(Bpp.cols(), Bpp.colidx(), Bpp.rowptr(), Bpp.values())?;

    // do P and Q iterations
    while !converged && i < max_it {
        // update iteration counter
        i = i + 1;

        // do P iteration, update Va //

        solver.solve(&LUp, &mut P, true)?;
        let dVa = neg(P.iter().copied());

        // update voltage
        for (&i, dVa) in zip(&pv_pq, dVa) {
            Va[i] += dVa;
        }
        V = polar(&Vm, &Va);

        // evaluate mismatch
        (P, Q) = {
            let Ibus: Vec<Complex64> = Ybus * &V;
            let Sbus = Sbus.s_bus(&Vm);
            // mis = (V .* conj(Ybus * V) - Sbus(Vm)) ./ Vm;
            let mis: Vec<Complex64> = izip!(&V, &Ibus, &Sbus, &Vm)
                .map(|(V, Ibus, Sbus, Vm)| V * Ibus.conj() - Sbus / Vm)
                .collect();
            (
                pv_pq.iter().map(|&i| mis[i].re).collect::<Vec<_>>(),
                pq.iter().map(|&i| mis[i].im).collect::<Vec<_>>(),
            )
        };

        // check tolerance
        let normP = norm_inf(&P);
        let normQ = norm_inf(&Q);
        log::debug!("normP = {}, normQ = {}", normP, normQ);
        if let Some(pm) = progress {
            pm.update(i, i - 1, normP, normQ);
        }
        if normP < tol && normQ < tol {
            converged = true;
            log::info!(
                "Fast-decoupled power flow converged in {} P-iterations and {} Q-iterations.",
                i,
                i - 1
            );
        }

        // do Q iteration, update Vm //

        solver.solve(&LUpp, &mut Q, true)?;
        let dVm = neg(Q.iter().copied());

        // update voltage
        zip(pq, dVm).for_each(|(&i, dVm)| Vm[i] += dVm);
        V = polar(&Vm, &Va);

        // evaluate mismatch
        (P, Q) = {
            let Ibus: Vec<Complex64> = Ybus * &V;
            let Sbus = Sbus.s_bus(&Vm);
            // mis = (V .* conj(Ybus * V) - Sbus(Vm)) ./ Vm;
            let mis: Vec<Complex64> = izip!(&V, &Ibus, &Sbus, &Vm)
                .map(|(V, Ibus, Sbus, Vm)| V * Ibus.conj() - Sbus / Vm)
                .collect();
            (
                pv_pq.iter().map(|&i| mis[i].re).collect::<Vec<_>>(),
                pq.iter().map(|&i| mis[i].im).collect::<Vec<_>>(),
            )
        };

        // check tolerance
        let normP = norm_inf(&P);
        let normQ = norm_inf(&Q);
        log::debug!("normP = {}, normQ = {}", normP, normQ);
        if let Some(pm) = progress {
            pm.update(i, i, normP, normQ);
        }
        if normP < tol && normQ < tol {
            converged = true;
            log::info!(
                "Fast-decoupled power flow converged in {} P-iterations and {} Q-iterations.",
                i,
                i
            );
        }
    }

    if !converged {
        log::info!(
            "Fast-decoupled power flow did not converge in {} iterations.",
            i
        );
    }

    Ok((V, converged, i))
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
) -> (CSR<usize, f64>, Option<CSR<usize, f64>>) {
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
        let (y_p, _) = make_ybus(base_mva, &bus, &branch, false);
        // y_p.map(|y| -y.im)
        -y_p.to_csr().imag()
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
        let (y_pp, _) = make_ybus(base_mva, &bus, &branch, false);
        // Some(y_pp.map(|y| -y.im))
        Some(-y_pp.to_csr().imag())
    } else {
        None
    };

    (b_p, b_pp)
}

pub trait ProgressMonitor {
    fn update(&self, p: usize, q: usize, norm_p: f64, norm_q: f64);
}

pub struct PrintProgress {}

impl ProgressMonitor for PrintProgress {
    fn update(&self, p: usize, q: usize, norm_p: f64, norm_q: f64) {
        if p == 0 {
            println!("iteration     max mismatch (p.u.)  ");
            println!("type   #        P            Q     ");
            println!("---- ----  -----------  -----------");
            // println!("  -  %3d   %10.3e   %10.3e", i, norm_p, norm_q);
            println!("  -  {}   {}   {}", p, norm_p, norm_q);
        } else {
            if p != q {
                // println!("  P  %3d   %10.3e   %10.3e", i, norm_p, norm_q);
                println!("  P  {}   {}   {}", p, norm_p, norm_q);
            } else {
                // println!("  Q  %3d   %10.3e   %10.3e", i, norm_p, norm_q);
                println!("  Q  {}   {}   {}", p, norm_p, norm_q);
            }
        }
    }
}
