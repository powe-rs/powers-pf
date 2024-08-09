use crate::dc::dc_pf;
use crate::fd;
use crate::gauss;
use crate::newton::*;
use crate::pfopt::{Alg, BusVoltage, GenQLimits, MPOpt, NodalBalance};
use crate::radial::radial_pf;

use powers::debug::format_polar_vec;
use powers::{bus_types, ext_to_int, int_to_ext, make_b_dc, make_sbus, make_ybus};
use powers::{MakeSBus, MPC};

use crate::pfmpc::PFMPC;
use crate::pfsoln::pfsoln;
use anyhow::Result;
use num_complex::Complex64;
use spsolve::{FactorSolver, Solver};
use std::collections::HashSet;
use std::f64::consts::PI;
use std::time::Instant;

pub fn runpf<F>(
    casedata: &MPC,
    mpopt: &MPOpt,
    solver: &dyn Solver<usize, f64>,
    factor_solver: &dyn FactorSolver<usize, f64, F>,
) -> Result<(MPC, bool)> {
    // options
    let qlim = mpopt.pf.enforce_q_limits != GenQLimits::IgnoreLimits; // enforce Q limits on gens?
    let dc = mpopt.dc; // use DC formulation?

    // read data
    let mpc = casedata;

    // convert to internal indexing
    let mut mpc = PFMPC::new(ext_to_int(mpc));
    let base_mva = mpc.base_mva;
    // let (base_mva, bus, gen, branch) = (mpc.base_mva, &mut mpc.bus, &mut mpc.gen, &mut mpc.branch);

    let (_t0, success, its) = if !mpc.bus.is_empty() {
        // get bus index lists of each type of bus
        let (ref_, pv, pq) = bus_types(&mpc.bus, &mpc.gen);

        //-----  run the power flow  -----
        let t0 = Instant::now();
        let mut success = false;
        let mut its = 0; // total iterations

        // if mpopt.verbose > 0 {
        //     v = mpver('all');
        //     fprintf('\nPowers Version % s, %s', v.Version, v.Date);
        // }
        if dc {
            // DC formulation
            // if mpopt.verbose > 0 {
            log::info!("DC Power Flow");
            // }
            // initial state
            let v_a0: Vec<f64> = mpc.bus.iter().map(|b| b.va * PI / 180.0).collect();

            // build B matrices and phase shift injections
            let (b_dc, b_f, p_businj, p_finj) = make_b_dc(base_mva, &mpc.bus, &mpc.branch);

            // compute complex bus power injections (generation - load)
            // adjusted for phase shifters and real shunts
            let s_bus = make_sbus(
                base_mva,
                &mpc.bus,
                &mpc.gen,
                mpopt.exp.sys_wide_zip_loads.pw,
                mpopt.exp.sys_wide_zip_loads.qw,
                None,
                None,
            );
            // let gs = bus.iter().map(|b| b.gs).collect_vec();
            // let Pbus: Vec<f64> = s_bus.real() - Pbusinj - gs / baseMVA;
            let p_bus: Vec<f64> = (0..mpc.bus.len())
                .map(|i| s_bus[i].re - p_businj[i] - mpc.bus[i].gs / base_mva)
                .collect();

            // "run" the power flow
            let (v_a, succ): (Vec<f64>, bool) =
                dc_pf(&b_dc, &p_bus, &v_a0, &ref_, &pv, &pq, solver)?;
            success = succ;

            its = 1;

            // update data matrices with solution
            let pf: Vec<f64> = (b_f * &v_a)
                .iter()
                .enumerate()
                .map(|(i, pf)| (pf + p_finj[i]) * base_mva)
                .collect();
            for (i, br) in mpc.branch.iter_mut().enumerate() {
                br.qf = Some(0.0);
                br.qt = Some(0.0);
                br.pf = Some(pf[i]);
                br.pt = Some(-pf[i]);
            }
            for (i, b) in mpc.bus.iter_mut().enumerate() {
                b.vm = 1.0;
                b.va = v_a[i] * 180.0 / PI;
            }
            // update Pg for slack generator (1st gen at ref bus)
            // (note: other gens at ref bus are accounted for in Pbus)
            //      Pg = Pinj + Pload + Gs
            //      newPg = oldPg + newPinj - oldPinj
            let b_ref = b_dc.select(Some(&ref_), None)?;
            let p_ref = b_ref * &v_a;
            for r in ref_ {
                for g in mpc.gen.iter_mut() {
                    if g.gen_bus == r {
                        g.pg += (p_ref[r] - p_bus[r]) * base_mva;
                        break;
                    }
                }
            }
        } else {
            let alg = mpopt.pf.algorithm;

            // initial state
            // let v0 = Arr::ones(bus.len()); // flat start
            let mut v0 = Vec::with_capacity(mpc.bus.len());
            for b in mpc.bus.iter() {
                v0.push(Complex64::from_polar(b.vm, b.va * PI / 180.0));
            }
            let pq_bus: HashSet<usize> = HashSet::from_iter(pq.clone().into_iter()); // exclude PQ buses
            for g in mpc.gen.iter() {
                if g.is_on() && !pq_bus.contains(&g.gen_bus) {
                    v0[g.gen_bus] =
                        Complex64::new(g.vg / (v0[g.gen_bus] * v0[g.gen_bus]).norm(), 0.0);
                }
            }
            log::debug!("V0: {}", format_polar_vec(&v0));

            if qlim {
                // let ref0 = ref_.clone(); // save index and angle of
                // let Varef0 = ref0.iter().map(|&r0| bus[r0].va).collect::<Vec<f64>>(); //   original reference bus(es)
                // let mut limited = vec![]; // list of indices of gens @ Q lims
                // let mut fixedQg = Arr::zeros(gen.len()); // Qg of gens at Q limits
            }

            // build admittance matrices
            let (y_bus, y_br) = make_ybus(base_mva, &mpc.bus, &mpc.branch, true);
            let (y_f, y_t) = y_br.unwrap();
            let (y_bus, y_f, y_t) = (y_bus.to_csr(), y_f.to_csr(), y_t.to_csr());
            log::trace!("Ybus:\n{}", y_bus.to_table());

            let mut repeat = true;
            while repeat {
                // function for computing V dependent complex bus power injections
                // (generation - load)
                // let s_bus: SBus = |v_m: &Arr<f64>| {
                //     make_sbus(baseMVA, &bus, &gen, mpopt, Some(v_m), None, false).0;
                // };
                let s_bus = MakeSBus {
                    base_mva,
                    bus: &mpc.bus,
                    gen: &mpc.gen,
                    pw: mpopt.exp.sys_wide_zip_loads.pw,
                    qw: mpopt.exp.sys_wide_zip_loads.qw,
                };

                let (v, succ, iterations) = match alg {
                    Alg::NR => {
                        let newtonpf_fcn = match mpopt.pf.current_balance {
                            NodalBalance::CURRENT => {
                                match mpopt.pf.v_cartesian {
                                    BusVoltage::POLAR => {
                                        newtonpf_i_polar // current, polar
                                    }
                                    BusVoltage::CARTESIAN => {
                                        newtonpf_i_cart // current, cartesian
                                    }
                                    BusVoltage::HYBRID => {
                                        newtonpf_i_hybrid // current, hybrid
                                    }
                                }
                            }
                            NodalBalance::POWER => {
                                match mpopt.pf.v_cartesian {
                                    BusVoltage::POLAR => {
                                        newtonpf_s_polar // default - power, polar
                                    }
                                    BusVoltage::CARTESIAN => {
                                        newtonpf_s_cart // power, cartesian
                                    }
                                    BusVoltage::HYBRID => {
                                        newtonpf_s_hybrid // power, hybrid
                                    }
                                }
                            }
                        };
                        newtonpf_fcn(&y_bus, &s_bus, &v0, &ref_, &pv, &pq, solver, &mpopt, None)?
                    }
                    Alg::FDBX | Alg::FDXB => {
                        let (b_p, b_pp) = fd::make_b(base_mva, &mpc.bus, &mpc.branch, alg, true);
                        let (b_p, b_pp) = (b_p.to_csr(), b_pp.unwrap().to_csr());
                        let progress = fd::PrintProgress {};
                        fd::fdpf(
                            &y_bus,
                            &s_bus,
                            &v0,
                            &b_p,
                            &b_pp,
                            &ref_,
                            &pv,
                            &pq,
                            factor_solver,
                            &mpopt,
                            Some(&progress),
                        )?
                    }
                    Alg::GS => {
                        let progress = gauss::PrintProgress {};
                        gauss::gausspf(
                            &y_bus,
                            &s_bus,
                            &v0,
                            &ref_,
                            &pv,
                            &pq,
                            solver,
                            &mpopt,
                            Some(&progress),
                        )?
                    }
                    Alg::SUM => {
                        let (_mpc, success, iterations) =
                            radial_pf(base_mva, &mpc.bus, &mpc.gen, &mpc.branch, mpopt)?;
                        (Vec::new(), success, iterations)
                    }
                };
                success = succ;
                its = its + iterations;

                // update data matrices with solution
                match alg {
                    Alg::NR | Alg::FDBX | Alg::FDXB | Alg::GS => {
                        pfsoln(
                            base_mva, &mut mpc,
                            // &mut mpc.bus,
                            // &mut mpc.gen,
                            // &mut mpc.branch,
                            &y_bus, &y_f, &y_t, &v, &ref_, &pv, &pq, &mpopt,
                        );
                        // bus = mpc_s.0;
                        // gen = mpc_s.1;
                        // branch = mpc_s.2;
                    }
                    _ => {}
                }

                if success && qlim {
                    // enforce generator Q limits
                    unimplemented!("generator Q limits");
                } else {
                    repeat = false; // don't enforce generator Q limits, once is enough
                }
            }
            // TODO: adjust voltage angles to make original ref bus correct
        }
        (t0, success, its)
    } else {
        // if mpopt.verbose {
        log::error!("Power flow not valid: Case contains no connected buses");
        // }
        (Instant::now(), false, 0)
    };
    // mpc.et = t0;
    mpc.success = Some(success);
    mpc.iterations = Some(its);

    // -----  output results  ----- //
    // convert back to original bus numbering & print results
    // mpc.bus = bus;
    // mpc.gen = gen;
    // mpc.branch = branch;
    let mut results = int_to_ext(&mpc).unwrap();

    let order = results.order().as_ref().unwrap().clone();

    // zero out result fields of out-of-service gens & branches
    let off = order.gen.status.off.clone();
    for i in off {
        results.gen[i].pg = 0.0;
        results.gen[i].qg = 0.0;
    }
    let off = order.branch.status.off.clone();
    for i in off {
        results.branch[i].pf = Some(0.0);
        results.branch[i].qf = Some(0.0);
        results.branch[i].pt = Some(0.0);
        results.branch[i].qt = Some(0.0);
    }

    // printpf(&results, 1, mpopt);

    Ok((results, success))
}

fn _have_zip_loads(mpopt: &MPOpt) -> bool {
    if let Some(pw) = mpopt.exp.sys_wide_zip_loads.pw {
        // TODO: check indexing
        if pw[1..2].iter().any(|&v| v != 0.0) {
            return true;
        }
    }
    if let Some(qw) = mpopt.exp.sys_wide_zip_loads.qw {
        if qw[1..2].iter().any(|&v| v != 0.0) {
            return true;
        }
    }
    false
    // (mpopt.exp.sys_wide_zip_loads.pw.is_some()
    //     && any(&mpopt.exp.sys_wide_zip_loads.pw.unwrap()[1..2]))
    //     || (mpopt.exp.sys_wide_zip_loads.qw.is_some()
    //         && any(&mpopt.exp.sys_wide_zip_loads.qw.unwrap()[1..2])) // TODO: check indexing
}

// fn printpf<W: Write>(_results: &MPC, _fd: W, _mpopt: &MPOpt) {}
