use crate::pfopt::MPOpt;
use caseformat::Gen;
use num_complex::Complex64;
use powers::total_load::{total_load, LoadType, LoadZone};
use powers::MPC;
use sparsetools::csr::CSR;
use std::collections::HashMap;
use std::f64::consts::PI;
use std::ops::Deref;

// Used in CPF.
pub fn pfsoln(
    base_mva: f64,
    mpc: &mut MPC,
    y_bus: &CSR<usize, Complex64>,
    y_f: &CSR<usize, Complex64>,
    y_t: &CSR<usize, Complex64>,
    v: &[Complex64],
    refbus: &[usize],
    _pv: &[usize],
    _pq: &[usize],
    mpopt: &MPOpt,
) {
    // let (mut bus, mut gen, mut branch) = (bus0.clone(), gen0.clone(), branch0.clone());
    // let (mut bus, mut gen, mut branch) = (Vec::from(bus0), Vec::from(gen0), Vec::from(branch0));
    // let mut bus: Vec<Bus> = Vec::from(bus0);

    for (i, b) in mpc.bus.iter_mut().enumerate() {
        b.vm = v[i].norm();
        b.va = v[i].arg() * 180.0 / PI;
    }

    let i_bus = y_bus * v;
    // let s_bus = gen.iter().map(|g| V[g.bus] * i_bus[g.bus].conj()).collect();

    let (_pd_gbus, qd_gbus) = total_load(
        &mpc.bus, /*(gbus, :)*/
        None,
        LoadZone::Bus,
        None,
        LoadType::Fixed,
        false,
        mpopt.exp.sys_wide_zip_loads.pw,
        mpopt.exp.sys_wide_zip_loads.qw,
        true,
    );
    let qd_gbus = qd_gbus.as_ref().unwrap();

    let is_on = |g: &Gen| g.is_on() && mpc.bus[g.gen_bus].is_pq();

    let nb = mpc.bus.len();
    let mut ngb = vec![0; nb];
    for g in mpc.gen.iter_mut() {
        let s_bus = v[g.gen_bus] * i_bus[g.gen_bus].conj(); // compute total injected bus power
        if is_on(&g.deref()) {
            g.qg = s_bus.im * base_mva + qd_gbus[g.gen_bus]; // inj Q + local Qd

            ngb[g.gen_bus] += 1;
        } else if g.is_off() {
            g.qg = 0.0;
        }
    }
    let ngg = mpc
        .gen
        .iter()
        .filter(|&g| is_on(g))
        .map(|g| ngb[g.gen_bus])
        .collect::<Vec<usize>>(); // number of gens at this gen's bus

    // ...at this point any buses with more than one generator will have
    // the total Q dispatch for the bus assigned to each generator. This
    // must be split between them. We do it first equally, then in proportion
    // to the reactive range of the generator.
    for (i, g) in mpc.gen.iter_mut().filter(|g| is_on(&g.deref())).enumerate() {
        g.qg = g.qg / ngg[i] as f64;

        let mut m = g.qg.abs();
        if !g.qmax.is_infinite() {
            m = m + g.qmax.abs();
        }
        if !g.qmin.is_infinite() {
            m = m + g.qmin.abs();
        }
        // TODO: each gen gets sum over all gens at same bus

        // replace +/- Inf limits with proxy +/- M
        if g.qmin.is_infinite() {
            g.qmin = if g.qmin.is_sign_positive() { m } else { -m };
        }
        if g.qmin.is_infinite() {
            g.qmax = if g.qmax.is_sign_positive() { m } else { -m };
        }
    }

    // Buses with more than one generator have the total reactive power
    // dispatch for the bus divided equally among each online generator.
    // The reactive power is divided in proportion to the reactive range
    // of each generator, according to the logic in the pfsoln function
    // by Ray Zimmerman from MATPOWER v7.
    let mut cg: HashMap<usize, Vec<usize>> = HashMap::new();
    for (i, g) in mpc.gen.iter().enumerate() {
        if is_on(g) {
            // cg[g.bus] = append(cg[g.bus], gen)
            if cg.contains_key(&g.gen_bus) {
                cg.get_mut(&g.gen_bus).unwrap().push(i)
            } else {
                cg.insert(g.gen_bus, vec![i]);
            }
            // match cg.get_mut(&g.bus) {
            //     Some(l) => l.push(g),
            //     None => cg.insert(g.bus, vec![g]),
            // }
        }
    }

    for l in cg.values() {
        if l.len() < 2 {
            continue;
        }
        let mut qg_tot = 0.0; // Total Qg at the bus.
        for &i in l.iter() {
            let g: &mut Gen = mpc.gen.get_mut(i).unwrap();
            qg_tot += g.qg;
        }

        // The sum of absolute Qg, Qmax and Qmin for all generators
        // at the bus (if Qmax/Qmin is non-infinite). Used as limit
        // when Qmax/Qmin is infinite.
        let mut m = 0.0;
        for &i in l.iter() {
            let g: &Gen = mpc.gen.get(i).unwrap();

            let mut mg = g.qg.abs();
            if !g.qmax.is_infinite() {
                mg = mg + g.qmax.abs();
            }
            if !g.qmin.is_infinite() {
                mg = mg + g.qmin.abs();
            }

            m += mg
        }
        let mut qmin = vec![0.0; l.len()];
        let mut qmax = vec![0.0; l.len()];

        for (i, &j) in l.iter().enumerate() {
            let g: &Gen = mpc.gen.get(j).unwrap();

            qmin[i] = g.qmin;
            qmax[i] = g.qmax;

            // replace +/- Inf limits with proxy +/- M
            if g.qmin.is_infinite() {
                qmin[i] = if g.qmin.is_sign_positive() { m } else { -m };
            }
            if g.qmax.is_infinite() {
                qmax[i] = if g.qmax.is_sign_positive() { m } else { -m };
            }
        }
        let qg_min: f64 = qmin.iter().sum(); // Minimum total Qg at the bus.
        let qg_max: f64 = qmax.iter().sum(); // Maximum total Qg at the bus.

        if (qg_min - qg_max).abs() > 1e-13 {
            let q = (qg_tot - qg_min) / (qg_max - qg_min);
            // for (i, mut pv) in l.iter_mut().enumerate() {
            for &i in l {
                let g: &mut Gen = mpc.gen.get_mut(i).unwrap();

                //pv.Qg = pv.QMin + (((QgTot - QgMin) / (QgMax - QgMin)) / (pv.QMax - pv.QMin))
                g.qg = g.qmin + (q * (g.qmax - g.qmin));
            }
        } else {
            // Zero Qg range at bus. Qg set such that all generators
            // at the bus violate their limits by the same amount.

            // Total mismatch at bus divided by number of online generators.
            let mis = (qg_tot - qg_min) / (l.len() as f64);
            for (i, &j) in l.iter().enumerate() {
                let g: &mut Gen = mpc.gen.get_mut(j).unwrap();

                g.qg = qmin[i] + mis;
            }
        }
    }

    // Update Pg for slack gen(s).
    for &r in refbus {
        let (pd_refk, _) = total_load(
            &vec![mpc.bus[r].clone()],
            None,
            LoadZone::Bus,
            None,
            LoadType::Fixed,
            false,
            mpopt.exp.sys_wide_zip_loads.pw,
            mpopt.exp.sys_wide_zip_loads.qw,
            false,
        );
        // let mut refgen = 0;
        let mut refgen0: Option<usize> = None;
        let mut sum = 0.0;
        for (i, g) in mpc.gen.iter_mut().enumerate() {
            if g.gen_bus == r {
                let s_bus = v[g.gen_bus] * i_bus[g.gen_bus].conj();

                g.pg = s_bus.re * base_mva + pd_refk[0]; // inj P + local Pd

                if refgen0.is_some() {
                    // more than one generator at this ref bus
                    sum += g.pg;
                }
                refgen0 = Some(i);
            }
        }
        // subtract off what is generated by other gens at this bus
        if let Some(i) = refgen0 {
            mpc.gen[i].pg -= sum;
        }
    }

    // Update/compute branch power flows.
    let i_fr_bus = y_f * v;
    let i_to_bus = y_t * v;
    for (i, br) in mpc.branch.iter_mut().enumerate() {
        if br.is_on() {
            let s_f: Complex64 = v[i] * i_fr_bus[i].conj() * base_mva; // complex power at "from" bus
            let s_t: Complex64 = v[i] * i_to_bus[i].conj() * base_mva; // complex power injected at "to" bus

            br.pf = Some(s_f.re);
            br.qf = Some(s_f.im);
            br.pt = Some(s_t.re);
            br.qt = Some(s_t.im);
        } else {
            br.pf = Some(0.0);
            br.qf = Some(0.0);
            br.pt = Some(0.0);
            br.qt = Some(0.0);
        }
    }

    // return (Vec::from(bus), Vec::from(gen), Vec::from(branch));
}
