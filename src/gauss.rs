use crate::pfopt::MPOpt;
use powers::SBus;

use anyhow::Result;
use full::slice::{abs, norm_inf};
use itertools::izip;
use num_complex::Complex64;
use powers::debug::format_rect_vec;
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

// Solves the power flow using a Gauss-Seidel method.
//
// Solves for bus voltages given the full system admittance matrix (for
// all buses), the complex bus power injection vector (for all buses),
// the initial vector of complex bus voltages, and column vectors with
// the lists of bus indices for the swing bus, PV buses, and PQ buses,
// respectively. The bus voltage vector contains the set point for
// generator (including ref bus) buses, and the reference angle of the
// swing bus, as well as an initial guess for remaining magnitudes and
// angles. MPOPT is a MATPOWER options struct which can be used to
// set the termination tolerance, maximum number of iterations, and
// output options (see MPOPTION for details). Uses default options
// if this parameter is not given. Returns the final complex voltages,
// a flag which indicates whether it converged or not, and the number
// of iterations performed.
#[allow(non_snake_case)]
pub(crate) fn gausspf(
    Ybus: &CSR<usize, Complex64>,
    SbusFn: &dyn SBus,
    V0: &[Complex64],
    _ref: &[usize],
    pv: &[usize],
    pq: &[usize],
    _solver: &dyn Solver<usize, f64>,
    mpopt: &MPOpt,
    progress: Option<&dyn ProgressMonitor>,
) -> Result<(Vec<Complex64>, bool, usize)> {
    let pv_pq = [pv, pq].concat();

    // options
    let tol = mpopt.pf.tolerance;
    let max_it = mpopt.pf.max_it_gs;

    // initialize
    let mut converged = false;
    let mut i = 0;
    let mut V = V0.to_vec();
    let Vm = abs(&V);

    // set up indexing for updating V
    let npv = pv.len();
    // let npq = pq.len();

    // evaluate F(x0)
    let Ibus: Vec<Complex64> = Ybus * &V;
    let mut Sbus = SbusFn.s_bus(&Vm);
    // let mis = &v * (y_bus * &v).conj() - s_bus.s_bus(&vm);
    // let mis: Vec<Complex64> = (0..nb).map(|i| v[i] * i_bus[i].conj() - s_bus[i]).collect();
    let mis: Vec<Complex64> = izip!(&V, &Ibus, &Sbus)
        .map(|(V, Ibus, Sbus)| V * Ibus.conj() - Sbus)
        .collect();
    let F: Vec<f64> = [
        pv_pq.iter().map(|&i| mis[i].re).collect::<Vec<_>>(),
        pq.iter().map(|&i| mis[i].im).collect::<Vec<_>>(),
    ]
    .concat();
    log::trace!("Sbus0: {}", format_rect_vec(&Sbus));

    // check tolerance
    let normF = norm_inf(&F);
    if let Some(pm) = progress {
        pm.update(i, normF);
    }
    if normF < tol {
        converged = true;
        log::info!("Converged!");
    }
    log::debug!("normF0: {}", normF);

    let diagYbus = Ybus.diagonal();

    // do Gauss-Seidel iterations
    while !converged && i < max_it {
        // update iteration counter
        i = i + 1;

        // update voltage
        // at PQ buses
        let Ibus: Vec<Complex64> = Ybus * &V;
        for &k in pq {
            // V(k) =  V(k) + (conj(Sbus(k) / V(k)) - Ybus(k,:) * V ) / Ybus(k,k);
            V[k] = V[k] + ((Sbus[k] / V[k]).conj() - Ibus[k]) / diagYbus[k];

            // + ((Sbus[k] / V[k]).conj() - Ybus.select(Some(&[k]), None)?.mat_vec(&V)?[0])
        }

        // at PV buses
        if npv != 0 {
            for &k in pv {
                // Sbus(k) = real(Sbus(k)) + 1j * imag( V(k) .* conj(Ybus(k,:) * V));
                Sbus[k] = Complex64::new(Sbus[k].re, (V[k] * Ibus[k].conj()).im);
                // V(k) =  V(k) + (conj(Sbus(k) / V(k)) - Ybus(k,:) * V ) / Ybus(k,k);
                V[k] = V[k] + ((Sbus[k] / V[k]).conj() - Ibus[k]) / diagYbus[k];

                // V(k) = Vm(k) * V(k) / abs(V(k));
            }
            for &k in pv {
                V[k] = Vm[k] * V[k] / Complex64::new(V[k].norm(), 0.0);
            }
        }

        // evalute F(x)
        let F = {
            let Ibus: Vec<Complex64> = Ybus * &V;
            let Sbus = SbusFn.s_bus(&Vm);
            log::trace!("Sbus_{}: {}", i, format_rect_vec(&Sbus));
            let mis: Vec<Complex64> = izip!(&V, &Ibus, &Sbus)
                .map(|(v, i_bus, s_bus)| v * i_bus.conj() - s_bus)
                .collect();
            [
                pv_pq.iter().map(|&i| mis[i].re).collect::<Vec<_>>(),
                pq.iter().map(|&i| mis[i].im).collect::<Vec<_>>(),
            ]
            .concat()
        };

        // check for convergence
        let norm_f = norm_inf(&F);
        if let Some(pm) = progress {
            pm.update(i, norm_f);
        }
        if norm_f < tol {
            converged = true;
            log::info!("Gauss-Seidel power flow converged in {} iterations.", i);
        }
        log::debug!("norm_f{}: {}", i, norm_f);
    }

    if !converged {
        log::info!(
            "Gauss-Seidel power flow did not converge in {} iterations.",
            i
        );
    }

    Ok((V, converged, i))
}
