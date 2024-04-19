use crate::newton::ProgressMonitor;
use crate::pfopt::MPOpt;
use powers::math::{norm_inf, J};
use powers::{d_imis_dv, SBus};

use anyhow::Result;
use num_complex::Complex64;
use sparsetools::coo::Coo;
use sparsetools::csr::{CCSR, CSR};
use spsolve::Solver;
use std::iter::zip;

/// Solves power flow using full Newton's method (current/cartesian).
///
/// Solves for bus voltages using a full Newton-Raphson method, using nodal
/// current balance equations and polar coordinate representation of
/// voltages.
///
/// The bus voltage vector contains the set point for generator
/// (including ref bus) buses, and the reference angle of the swing
/// bus, as well as an initial guess for remaining magnitudes and
/// angles.
///
/// Returns the final complex voltages, a flag which indicates whether it
/// converged or not, and the number of iterations performed.
pub(crate) fn newtonpf_i_polar(
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
    let pv_pq = [pv, pq].concat();

    let tol = mpopt.pf.tolerance;
    let max_it = mpopt.pf.max_it_nr;

    let mut converged = false;
    let mut i = 0;
    let mut v = v0.to_vec();
    let mut va: Vec<f64> = v.iter().map(|v| v.arg()).collect();
    let mut vm: Vec<f64> = v.iter().map(|v| v.norm()).collect();
    // let vm_pv: Vec<f64> = pv.iter().map(|&i| vm[i]).collect();
    let nb = v0.len();

    // set up indexing for updating V
    let npv = pv.len();
    let npq = pq.len();
    let (j1, j2) = (0, npv); // j1:j2 - V angle of pv buses
    let (j3, j4) = (j2, j2 + npq); // j3:j4 - V angle of pq buses
    let (j5, j6) = (j4, j4 + npv); // j5:j6 - Q of pv buses
    let (j7, j8) = (j6, j6 + npv); // j7:j8 - V mag of pq buses
    let j1_j2: Vec<usize> = (j1..j2).collect();
    let j3_j4: Vec<usize> = (j3..j4).collect();
    let j5_j6: Vec<usize> = (j5..j6).collect();
    let j7_j8: Vec<usize> = (j7..j8).collect();

    // evaluate F(x0)
    let mut s_bus = s_bus_fn.s_bus(&vm);
    let mis: Vec<Complex64> = {
        // let ybus_pv = y_bus.select(Some(&pv), None)?;
        // let i_bus_pv: Vec<Complex64> = ybus_pv * &v;
        let i_bus: Vec<Complex64> = y_bus * &v;
        // for i in 0..npv {
        pv.iter().for_each(|&i| {
            s_bus[i].im = (v[i] * i_bus[i]).conj().im;
        });
        // (0..nb)
        //     .map(|i| i_bus[i] - (s_bus[i] / v[i]).conj())
        //     .collect()
        zip(&v, zip(&i_bus, &s_bus))
            .map(|(i_bus, (s_bus, v))| i_bus - (s_bus / v).conj())
            .collect()
    };
    let f = [
        pv_pq.iter().map(|&i| mis[i].re).collect::<Vec<_>>(),
        pv_pq.iter().map(|&i| mis[i].im).collect::<Vec<_>>(),
    ]
    .concat();

    // check tolerance
    let norm_f = norm_inf(&f).unwrap();
    if let Some(pm) = progress {
        pm.update(i, norm_f); // max Ir & Ii mismatch (p.u.)
    }
    if norm_f < tol {
        converged = true;
        log::info!("Converged!");
    }

    fn compose(j: [[&Coo<usize, f64>; 3]; 2]) -> Result<Coo<usize, f64>> {
        let j1x = Coo::h_stack3(j[0][0], j[0][1], j[0][2])?;
        let j2x = Coo::h_stack3(j[1][0], j[1][1], j[1][2])?;
        Coo::v_stack(&j1x, &j2x)
    }

    // do Newton iterations
    while !converged && i < max_it {
        // update iteration counter
        i = i + 1;

        // evaluate Jacobian
        /*
        dImis_dQ = sparse(pv, pv, 1j./conj(V(pv)), n, n);
        [dImis_dVa, dImis_dVm] = dImis_dV(Sb, Ybus, V);
        dImis_dVm(:, pv) = dImis_dQ(:, pv);

        j11 = real(dImis_dVa([pv; pq], [pv; pq]));
        j12 = real(dImis_dVm([pv; pq], [pv; pq]));
        j21 = imag(dImis_dVa([pv; pq], [pv; pq]));
        j22 = imag(dImis_dVm([pv; pq], [pv; pq]));

        J = [   j11 j12;
                j21 j22;    ];
        */
        let d_imis_d_q = Coo::new(
            nb,
            nb,
            pv.to_vec(),
            pv.to_vec(),
            pv.iter().map(|&i| J / v[i].conj()).collect(),
        )?
        .to_csr();
        let (d_imis_d_va, d_imis_d_vm) = d_imis_dv(&s_bus, &y_bus, &v, false)?;
        // TODO: dImis_dVm(:, pv) = dImis_dQ(:, pv);

        /*
        j11 = real(dImis_dVa([pv; pq], [pv; pq]));
        j12 = real(dImis_dQ([pv; pq], pv));
        j13 = real(dImis_dVm([pv; pq], pq));
        j21 = imag(dImis_dVa([pv; pq], [pv; pq]));
        j22 = imag(dImis_dQ([pv; pq], pv));
        j23 = imag(dImis_dVm([pv; pq], pq));

        J = [   j11 j12 j13;
                j21 j22 j23;    ];
        */
        let j11 = d_imis_d_va.select(Some(&pv_pq), Some(&pv_pq))?.real();
        let j12 = d_imis_d_q.select(Some(&pv_pq), Some(pv))?.real();
        let j13 = d_imis_d_vm.select(Some(&pv_pq), Some(pq))?.real();
        let j21 = d_imis_d_va.select(Some(&pv_pq), Some(&pv_pq))?.imag();
        let j22 = d_imis_d_q.select(Some(&pv_pq), Some(pv))?.imag();
        let j23 = d_imis_d_vm.select(Some(&pv_pq), Some(pq))?.imag();

        let jac = compose([
            [&j11.to_coo(), &j12.to_coo(), &j13.to_coo()],
            [&j21.to_coo(), &j22.to_coo(), &j23.to_coo()],
        ])?
        .to_csc();

        // compute update step
        // let dx = lin_solver.solve(jac.view(), f.iter().map(|&f_i| -f_i).collect())?;
        let dx = {
            let mut neg_f: Vec<f64> = f.iter().map(|f| -f).collect();
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

        // for (i, j) in (j1..j2).enumerate() {
        //     va[pv[i]] += dx[j];
        // }
        // for (i, j) in (j3..j4).enumerate() {
        //     va[pq[i]] += dx[j];
        // }
        // for (i, j) in (j5..j6).enumerate() {
        //     vm[pq[i]] += dx[j];
        // }

        // update voltage
        if npv != 0 {
            /*
            Va(pv) = Va(pv) + dx(j1:j2);
            Sb(pv) = real(Sb(pv)) + 1j * (imag(Sb(pv)) + dx(j5:j6));
            */
            // let dx_va = dx.select(&(j1..j2).collect::<Vec<usize>>());
            // let dx_sb = dx.select(&(j5..j6).collect::<Vec<usize>>());
            // for (i, &j) in pv.iter().enumerate() {
            //     va[j] += dx_va[i];
            //     s_bus[j] += Complex64::new(0.0, dx_sb[i]);
            // }

            (0..npv).for_each(|i| va[pv[i]] += dx[j1_j2[i]]);
            (0..npv).for_each(|i| s_bus[pv[i]].im += dx[j5_j6[i]]);
            // for i in 0..npv {
            //     va[pv[i]] += dx[j1_j2[i]];
            //     s_bus[pv[i]].im += dx[j5_j6[i]];
            // }
        }
        if npq != 0 {
            /*
            Va(pq) = Va(pq) + dx(j3:j4);
            Vm(pq) = Vm(pq) + dx(j7:j8);
            */
            (0..npq).for_each(|i| va[pq[i]] += dx[j3_j4[i]]);
            (0..npq).for_each(|i| vm[pq[i]] += dx[j7_j8[i]]);
            // for i in 0..npq {
            //     va[pq[i]] += dx[j3_j4[i]];
            //     vm[pq[i]] += dx[j7_j8[i]];
            // }
        }

        // update Vm and Va again in case we wrapped around with a negative Vm
        // v = Vec::from_polar(&vm, &va);
        v = zip(vm, va)
            .map(|(vm, va)| Complex64::from_polar(vm, va))
            .collect();
        va = v.iter().map(|v| v.arg()).collect();
        vm = v.iter().map(|v| v.norm()).collect();

        // evalute F(x)
        // mis = Ybus * V - conj(Sb ./ V);
        let i_bus = y_bus * &v;
        let mis: Vec<_> = zip(&v, zip(&i_bus, &s_bus))
            .map(|(i_bus, (s_bus, v))| i_bus - (s_bus / v).conj())
            .collect();
        // let mis: Vec<Complex64> = (0..nb)
        //     .map(|i| i_bus[i] - (s_bus[i] / v[i]).conj())
        //     .collect();

        let f = [
            pv_pq.iter().map(|&i| mis[i].re).collect::<Vec<f64>>(),
            pv_pq.iter().map(|&i| mis[i].im).collect::<Vec<f64>>(),
        ]
        .concat();

        // check for convergence
        let norm_f = norm_inf(&f).unwrap();
        if let Some(pm) = progress {
            pm.update(i, norm_f);
        }
        if norm_f < tol {
            converged = true;
            log::info!(
                "Newton's method power flow (current balance, polar) converged in {} iterations.",
                i
            );
        }
    }

    if !converged {
        log::info!(
            "Newton's method power flow (current balance, polar) did not converge in {} iterations.",
            i
        );
    }

    Ok((v, converged, i))
}

/// Solves power flow using full Newton's method (current/cartesian).
///
/// Solves for bus voltages using a full Newton-Raphson method, using nodal
/// current balance equations and cartesian coordinate representation of
/// voltages.
///
/// The bus voltage vector contains the set point for generator
/// (including ref bus) buses, and the reference angle of the swing
/// bus, as well as an initial guess for remaining magnitudes and
/// angles.
///
/// Returns the final complex voltages, a flag which indicates whether it
/// converged or not, and the number of iterations performed.
pub(crate) fn newtonpf_i_cart(
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
    let pv_pq = [pv, pq].concat();

    let tol = mpopt.pf.tolerance;
    let max_it = mpopt.pf.max_it_nr;

    let mut converged = false;
    let mut i = 0;
    let mut v = v0.to_vec();
    // let mut va: Vec<f64> = v.iter().map(|v| v.arg()).collect();
    let vm: Vec<f64> = v.iter().map(|v| v.norm()).collect();
    // let vm_pv = pv.iter().map(|&i| vm[i]).collect();
    let nb = v0.len();

    // set up indexing for updating V
    let npv = pv.len();
    let npq = pq.len();
    let (j1, j2) = (0, npv); // j1:j2 - Q of pv buses
    let (j3, j4) = (j2, j2 + npq); // j3:j4 - Vr of pq buses
    let (j5, j6) = (j4, j4 + npv); // j5:j6 - Vr of pv buses
    let (j7, j8) = (j6, j6 + npq); // j7:j9 - Vi of pq buses
    let (j9, j10) = (j8, j8 + npv); // j9:j10 - Vi of pv buses
    let j1_j2: Vec<usize> = (j1..j2).collect();
    let j3_j4: Vec<usize> = (j3..j4).collect();
    let j5_j6: Vec<usize> = (j5..j6).collect();
    let j7_j8: Vec<usize> = (j7..j8).collect();
    let j9_j10: Vec<usize> = (j9..j10).collect();

    // evaluate F(x0)
    let mut s_bus = s_bus_fn.s_bus(&vm);
    let mis: Vec<_> = {
        // let ybus_pv = y_bus.select(Some(&pv), None)?;
        // let i_bus_pv: Vec<Complex64> = ybus_pv * &v;
        // for i in 0..npv {
        //     s_bus[pv[i]].im = (v[pv[i]] * i_bus_pv[i]).conj().im;
        // }
        let i_bus: Vec<Complex64> = y_bus * &v;
        pv.iter()
            .for_each(|&i| s_bus[i].im = (v[i] * i_bus[i]).conj().im);

        zip(&v, zip(&i_bus, &s_bus))
            .map(|(i_bus, (s_bus, v))| i_bus - (s_bus / v).conj())
            .collect()
    };

    let f: Vec<_> = [
        pv_pq.iter().map(|&i| mis[i].re).collect::<Vec<_>>(),
        pv_pq.iter().map(|&i| mis[i].im).collect(),
        pv.iter()
            .map(|&i| (v[i] * v[i].conj()).re - (vm[i] * vm[i]))
            .collect(),
    ]
    .concat();

    // let sb_pv = Vec::from_parts(
    //     &s_bus.select(&pv).real(),
    //     &(v.select(&pv) * conj(&(y_bus * &v))).imag(),
    // );
    // s_bus.set(pv, &sb_pv);
    // let mis = (y_bus * &v) - conj(&(&s_bus / &v));
    //
    //     let mut f = {
    //         let mis_pvpq = mis.select(&pv_pq);
    //         let v_pv = v.select(&pv);
    //         Vec::concat(&[
    //             &mis_pvpq.real(),
    //             &mis_pvpq.imag(),
    //             &(&(&v_pv * conj(&v_pv)).real() - vm_pv.pow(2)),
    //         ])
    //     };

    // check tolerance
    let norm_f = norm_inf(&f).unwrap();
    if let Some(pm) = progress {
        pm.update(i, norm_f); // max Ir & Ii mismatch (p.u.)
    }
    if norm_f < tol {
        converged = true;
        log::info!("Converged!");
    }

    // do Newton iterations
    while !converged && i < max_it {
        // update iteration counter
        i = i + 1;

        // evaluate Jacobian
        /*
        dImis_dQ = sparse(pv, pv, 1j./conj(V(pv)), n, n);
        dV2_dVr = sparse(1:npv, npq+(1:npv), 2*real(V(pv)), npv, npv+npq);
        dV2_dVi = sparse(1:npv, npq+(1:npv), 2*imag(V(pv)), npv, npv+npq);
        [dImis_dVr, dImis_dVi] = dImis_dV(Sb, Ybus, V, 1);
             */
        // let v_pv = v.select(&pv);
        let d_imis_d_q = Coo::new(
            nb,
            nb,
            pv.to_vec(),
            pv.to_vec(),
            pv.iter().map(|&i| J / v[i].conj()).collect(),
        )?
        .to_csr();
        let d_v2_d_vr = Coo::new(
            npv,
            npv + npq,
            (0..npv).collect(),
            (npq..npq + npv).collect(),
            pv.iter().map(|&i| 2.0 * v[i].re).collect(),
        )?
        .to_csr();
        let d_v2_d_vi = Coo::new(
            npv,
            npv + npq,
            (0..npv).collect(),
            (npq..npq + npv).collect(),
            pv.iter().map(|&i| 2.0 * v[i].im).collect(),
        )?
        .to_csr();
        let (d_imis_d_vr, d_imis_d_vi) = d_imis_dv(&s_bus, &y_bus, &v, true)?;

        // handling of derivatives for voltage dependent loads
        // (not yet implemented) goes here

        let j11 = d_imis_d_q.select(Some(pv), Some(pv))?.real();
        let j12 = d_imis_d_vr.select(Some(&pv_pq), Some(&pv_pq))?.real();
        let j13 = d_imis_d_vi.select(Some(&pv_pq), Some(&pv_pq))?.real();

        let j21 = d_imis_d_q.select(Some(&pv_pq), Some(&pv))?.imag();
        let j22 = d_imis_d_vr.select(Some(&pv_pq), Some(&pv_pq))?.imag();
        let j23 = d_imis_d_vi.select(Some(&pv_pq), Some(&pv_pq))?.imag();

        let j31 = CSR::with_size(npv, npv);
        let j32 = d_v2_d_vr;
        let j33 = d_v2_d_vi;
        /*
        j11 = real(dImis_dQ([pv; pq], pv));
        j12 = real(dImis_dVr([pv; pq], [pq; pv]));
        j13 = real(dImis_dVi([pv; pq], [pq; pv]));

        j21 = imag(dImis_dQ([pv; pq], pv));
        j22 = imag(dImis_dVr([pv; pq], [pq; pv]));
        j23 = imag(dImis_dVi([pv; pq], [pq; pv]));

        j31 = sparse(npv, npv);
        j32 = dV2_dVr;
        j33 = dV2_dVi;
        */

        let jac = Coo::compose3([
            [&j11.to_coo(), &j12.to_coo(), &j13.to_coo()],
            [&j21.to_coo(), &j22.to_coo(), &j23.to_coo()],
            [&j31.to_coo(), &j32.to_coo(), &j33.to_coo()],
        ])?
        .to_csc();

        // compute update step
        let dx = {
            let mut neg_f: Vec<f64> = f.iter().map(|f| -f).collect();
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

        // update voltage
        if npv != 0 {
            /*
            V(pv) = V(pv) + dx(j5:j6) + 1j * dx(j9:j10);
            Sb(pv) = real(Sb(pv)) + 1j * (imag(Sb(pv)) + dx(j1:j2));
            */
            // let dx_v: Vec<Complex64> = Vec::from_parts(
            //     &dx.select(&(j5..j6).collect::<Vec<usize>>()),
            //     &dx.select(&(j9..j10).collect::<Vec<usize>>()),
            // );
            // let dx_sb = dx.select(&(j1..j2).collect::<Vec<usize>>());
            // for (i, &j) in pv.iter().enumerate() {
            //     v[j] += dx_v[i];
            //     s_bus[j] += Complex64::new(0.0, dx_sb[i]);
            // }

            for i in 0..npv {
                v[pv[i]] += Complex64::new(dx[j5_j6[i]], dx[j9_j10[i]]);
                s_bus[pv[i]].im += dx[j1_j2[i]];
            }
        }
        if npq != 0 {
            // V(pq) = V(pq) + dx(j3:j4) + 1j * dx(j7:j8);
            // let dx_v = Vec::<Complex64>::from_parts(
            //     &dx.select(&(j3..j4).collect::<Vec<usize>>()),
            //     &dx.select(&(j7..j8).collect::<Vec<usize>>()),
            // );
            // for (i, &j) in pq.iter().enumerate() {
            //     v[j] += dx_v[i];
            // }

            for i in 0..npq {
                v[pq[i]] += Complex64::new(dx[j3_j4[i]], dx[j7_j8[i]]);
            }
        }

        // evalute F(x)
        /*
        mis = Ybus * V - conj(Sb ./ V);
        F = [   real(mis([pv; pq]));
                imag(mis([pv; pq]));
                V(pv) .* conj(V(pv)) - Vmpv.^2  ];
        */
        let i_bus = y_bus * &v;
        // let mis = (0..nb)
        //     .map(|i| i_bus[i] - (s_bus[i] / v[i]).conj())
        //     .collect::<Vec<Complex64>>();
        let mis: Vec<_> = zip(&v, zip(&i_bus, &s_bus))
            .map(|(i_bus, (s_bus, v))| i_bus - (s_bus / v).conj())
            .collect();
        let f: Vec<_> = [
            pv_pq.iter().map(|&i| mis[i].re).collect::<Vec<_>>(),
            pv_pq.iter().map(|&i| mis[i].im).collect(),
            pv.iter()
                .map(|&i| (v[i] * v[i].conj()).re - (vm[i] * vm[i]))
                .collect(),
        ]
        .concat();

        // check for convergence
        let norm_f = norm_inf(&f).unwrap();
        if let Some(pm) = progress {
            pm.update(i, norm_f);
        }
        if norm_f < tol {
            converged = true;
            log::info!(
                "Newton's method power flow (current balance, cartesian) converged in {} iterations.",
                i
            );
        }
    }

    if !converged {
        log::info!(
            "Newton's method power flow (current balance, cartesian) did not converge in {} iterations.",
            i
        );
    }

    Ok((v, converged, i))
}

/// Solves power flow using full Newton's method (current/hybrid).
///
/// Solves for bus voltages using a full Newton-Raphson method, using nodal
/// current balance equations and a hybrid representation of voltages, where
/// a polar update is computed using a cartesian Jacobian.
///
/// The bus voltage vector contains the set point for generator
/// (including ref bus) buses, and the reference angle of the swing
/// bus, as well as an initial guess for remaining magnitudes and
/// angles.
///
/// Returns the final complex voltages, a flag which indicates whether it
/// converged or not, and the number of iterations performed.
pub(crate) fn newtonpf_i_hybrid(
    y_bus: &CSR<usize, Complex64>,
    s_bus_fn: &dyn SBus,
    v0: &[Complex64],
    _ref: &[usize],
    pv: &[usize],
    pq: &[usize],
    _solver: &dyn Solver<usize, f64>,
    mpopt: &MPOpt,
    progress: Option<&dyn ProgressMonitor>,
) -> Result<(Vec<Complex64>, bool, usize)> {
    let pv_pq = [pv, pq].concat();

    let tol = mpopt.pf.tolerance;
    let max_it = mpopt.pf.max_it_nr;

    let mut converged = false;
    let mut i = 0;
    let v = v0.to_vec();
    // let mut va = v.iter().map(|v| v.arg()).collect_vec();
    let vm: Vec<_> = v.iter().map(|v| v.norm()).collect();
    // let vm_pv = vm.select(pv);
    // let n = v0.len();

    // set up indexing for updating V
    let npv = pv.len();
    let npq = pq.len();
    let (_j1, j2) = (0, npv); // j1:j2 - Q of pv buses
    let (_j3, j4) = (j2, j2 + npq); // j3:j4 - Vr of pq buses
    let (_j5, j6) = (j4, j4 + npv); // j5:j6 - Vi of pv buses
    let (_j7, _j8) = (j6, j6 + npq); // j7:j8 - Vi of pq buses

    // evaluate F(x0)
    let mut s_bus = s_bus_fn.s_bus(&vm);
    // Sb(pv) = real(Sb(pv)) + 1j * imag(V(pv) .* conj(Ybus(pv, :) * V));
    let mis: Vec<_> = {
        // let ybus_pv = y_bus.select(Some(&pv), None)?;
        // let i_bus_pv: Vec<Complex64> = ybus_pv * &v;
        // for i in 0..npv {
        //     s_bus[pv[i]].im = (v[pv[i]] * i_bus_pv[i]).conj().im;
        // }
        let i_bus: Vec<Complex64> = y_bus * &v;
        pv.iter()
            .for_each(|&i| s_bus[i].im = (v[i] * i_bus[i]).conj().im);

        zip(&v, zip(&i_bus, &s_bus))
            .map(|(i_bus, (s_bus, v))| i_bus - (s_bus / v).conj())
            .collect()
        // (0..nb)
        //     .map(|i| i_bus[i] - (s_bus[i] / v[i]).conj())
        //     .collect()
    };

    let f = [
        pv_pq.iter().map(|&i| mis[i].re).collect::<Vec<_>>(),
        pv_pq.iter().map(|&i| mis[i].im).collect::<Vec<_>>(),
    ]
    .concat();

    // check tolerance
    let norm_f = norm_inf(&f).unwrap();
    if let Some(pm) = progress {
        pm.update(i, norm_f); // max Ir & Ii mismatch (p.u.)
    }
    if norm_f < tol {
        converged = true;
        log::info!("Converged!");
    }

    // do Newton iterations
    while !converged && i < max_it {
        // update iteration counter
        i = i + 1;

        // evaluate Jacobian
        /*
        dImis_dQ = sparse(1:n, 1:n, 1j./conj(V), n, n);
        [dImis_dVr, dImis_dVi] = dImis_dV(Sb, Ybus, V, 1);
        dImis_dVi(:, pv) = ...
            dImis_dVi(:, pv) * sparse(1:npv, 1:npv, real(V(pv)), npv, npv) - ...
            dImis_dVr(:, pv) * sparse(1:npv, 1:npv, imag(V(pv)), npv, npv);
        dImis_dVr(:, pv) = dImis_dQ(:, pv);
        */

        // let d_imis_d_q = Coo::with_diagonal(v.iter().map(|&vi| J / vi.conj()).collect()).to_csr();
        // let (d_imis_d_vr, d_imis_d_vi) = d_imis_d_v(&s_bus, &y_bus, &v, true)?;
    }

    unimplemented!("newtonpf_i_hybrid")
}
