use crate::newton::ProgressMonitor;
use crate::pfopt::MPOpt;
use powers::debug::{format_f64_vec, format_polar_vec, format_rect_vec};
use powers::{d_sbus_d_v, SBus};

use anyhow::Result;
use full::slice::norm_inf;
use num_complex::Complex64;
use sparsetools::coo::Coo;
use sparsetools::csr::{CCSR, CSR};
use spsolve::Solver;
use std::iter::zip;

/// Solves power flow using full Newton's method (power/polar).
///
/// Solves for bus voltages using a full Newton-Raphson method, using nodal
/// power balance equations and polar coordinate representation of
/// voltages.
///
/// The bus voltage vector contains the set point for generator
/// (including ref bus) buses, and the reference angle of the swing
/// bus, as well as an initial guess for remaining magnitudes and
/// angles.
///
/// Returns the final complex voltages, a flag which indicates whether it
/// converged or not, and the number of iterations performed.
pub(crate) fn newtonpf_s_polar(
    y_bus: &CSR<usize, Complex64>,
    s_bus_fn: &dyn SBus,
    v0: &[Complex64],
    _ref: &[usize],
    pv: &[usize],
    pq: &[usize],
    solver: &dyn Solver<usize, f64>,
    mpopt: &MPOpt,
    progress: Option<&dyn ProgressMonitor>,
) -> Result<(Vec<Complex64>, bool, usize)> {
    // let nb = v0.len();
    let pv_pq = [pv, pq].concat();

    let tol = mpopt.pf.tolerance;
    let max_it = mpopt.pf.max_it_nr;
    // let lin_solver  = mpopt.pf.nr.lin_solver;

    let mut converged = false;
    let mut i = 0;
    let mut v: Vec<Complex64> = v0.to_vec();
    let mut va: Vec<f64> = v.iter().map(|v| v.arg()).collect();
    let mut vm: Vec<f64> = v.iter().map(|v| v.norm()).collect();

    // set up indexing for updating V
    let npv = pv.len();
    let npq = pq.len();
    let (j1, j2) = (0, npv); // j1:j2 - V angle of pv buses
    let (j3, j4) = (j2, j2 + npq); // j3:j4 - V angle of pq buses
    let (j5, j6) = (j4, j4 + npq); // j5:j6 - V mag of pq buses

    // let j1_j2 = j1..j2;
    // let j3_j4 = j3..j4;
    // let j5_j6 = j5..j6;

    // let j1_j2 = 0..npv; // V angle of pv buses
    // let j3_j4 = j1_j2.end..(j1_j2.end + npq); // V angle of pq buses
    // let j5_j6 = j3_j4.end..(j3_j4.end + npq); // V mag of pq buses

    // evaluate F(x0)
    let i_bus: Vec<Complex64> = y_bus * &v;
    let s_bus = s_bus_fn.s_bus(&vm);
    // let mis = &v * (y_bus * &v).conj() - s_bus.s_bus(&vm);
    // let mis: Vec<Complex64> = (0..nb).map(|i| v[i] * i_bus[i].conj() - s_bus[i]).collect();
    let mis: Vec<Complex64> = zip(&v, zip(&i_bus, &s_bus))
        .map(|(v, (i_bus, s_bus))| v * i_bus.conj() - s_bus)
        .collect();
    let mut f: Vec<f64> = [
        pv_pq.iter().map(|&i| mis[i].re).collect::<Vec<_>>(),
        pq.iter().map(|&i| mis[i].im).collect::<Vec<_>>(),
    ]
    .concat();
    log::trace!("Sbus0: {}", format_rect_vec(&s_bus));

    // check tolerance
    let norm_f = norm_inf(&f);
    if let Some(pm) = progress {
        pm.update(i, norm_f);
    }
    if norm_f < tol {
        converged = true;
        log::info!("Converged!");
    }
    log::debug!("norm_f0: {}", norm_f);

    // do Newton iterations
    while !converged && i < max_it {
        // update iteration counter
        i = i + 1;

        // evaluate Jacobian
        let (d_sbus_d_va, d_sbus_d_vm) = d_sbus_d_v(y_bus, &v, false);
        log::trace!("dSbus_dVa:\n{}", d_sbus_d_va.to_table());
        log::trace!("dSbus_dVm:\n{}", d_sbus_d_vm.to_table());
        let neg_d_sd_d_vm = s_bus_fn.d_sbus_d_vm(&vm);
        let d_sbus_d_vm = d_sbus_d_vm - neg_d_sd_d_vm;

        let j11 = d_sbus_d_va.select(Some(&pv_pq), Some(&pv_pq))?.real();
        let j12 = d_sbus_d_vm.select(Some(&pv_pq), Some(pq))?.real();
        let j21 = d_sbus_d_va.select(Some(pq), Some(&pv_pq))?.imag();
        let j22 = d_sbus_d_vm.select(Some(pq), Some(pq))?.imag();

        let jac = Coo::compose([
            [&j11.to_coo(), &j12.to_coo()],
            [&j21.to_coo(), &j22.to_coo()],
        ])?
        .to_csc();
        log::trace!("J_{}:\n{}", i, jac.to_csr().to_table());

        // compute update step
        let dx = {
            let mut neg_f: Vec<f64> = f.iter().map(|f| -f).collect();
            log::trace!("-F: {}", format_f64_vec(&neg_f));
            solver.solve(
                jac.cols(),
                jac.rowidx(),
                jac.colptr(),
                jac.values(),
                &mut neg_f,
                false,
            )?;
            neg_f
        };
        log::trace!("dx: {}", format_f64_vec(&dx));

        // update voltage
        // pv.iter().zip(j1..j2).for_each(|(&i, j)| va[i] += dx[j]);
        // pq.iter().zip(j3..j4).for_each(|(&i, j)| va[i] += dx[j]);
        // pq.iter().zip(j5..j6).for_each(|(&i, j)| vm[i] += dx[j]);
        for (i, j) in (j1..j2).enumerate() {
            va[pv[i]] += dx[j];
        }
        for (i, j) in (j3..j4).enumerate() {
            va[pq[i]] += dx[j];
        }
        for (i, j) in (j5..j6).enumerate() {
            vm[pq[i]] += dx[j];
        }

        // update Vm and Va again in case we wrapped around with a negative Vm
        v = zip(vm, va)
            .map(|(vm, va)| Complex64::from_polar(vm, va))
            .collect();
        va = v.iter().map(|v| v.arg()).collect();
        vm = v.iter().map(|v| v.norm()).collect();
        log::debug!("V_{}: {}", i, format_polar_vec(&v));

        // evalute F(x)
        let i_bus: Vec<Complex64> = y_bus * &v;
        let s_bus = s_bus_fn.s_bus(&vm);
        // let mis = &v * conj(&(y_bus * &v)) - s_bus.s_bus(&vm);
        // let mis: Vec<Complex64> = (0..nb).map(|i| v[i] * i_bus[i].conj() - s_bus[i]).collect();
        let mis: Vec<Complex64> = zip(&v, zip(&i_bus, &s_bus))
            .map(|(v, (i_bus, s_bus))| v * i_bus.conj() - s_bus)
            .collect();
        f = [
            pv_pq.iter().map(|&i| mis[i].re).collect::<Vec<_>>(),
            pq.iter().map(|&i| mis[i].im).collect::<Vec<_>>(),
        ]
        .concat();
        log::trace!("Sbus_{}: {}", i, format_rect_vec(&s_bus));

        // check for convergence
        let norm_f = norm_inf(&f);
        if let Some(pm) = progress {
            pm.update(i, norm_f);
        }
        if norm_f < tol {
            converged = true;
            log::info!(
                "Newton's method power flow (power balance, polar) converged in {} iterations.",
                i
            );
        }
        log::debug!("norm_f{}: {}", i, norm_f);
    }

    if !converged {
        log::info!(
            "Newton's method power flow (power balance, polar) did not converge in {} iterations.",
            i
        );
    }

    Ok((v, converged, i))
}

pub(crate) fn newtonpf_s_cart(
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
    unimplemented!("newtonpf_s_cart");
}

pub(crate) fn newtonpf_s_hybrid(
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
    unimplemented!("newtonpf_s_hybrid")
}
