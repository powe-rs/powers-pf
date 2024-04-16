use anyhow::Result;
use powers::math::norm_inf;
use sparsetools::csr::CSR;
use spsolve::Solver;
use std::iter::zip;

/// Solves a DC power flow.
///
/// Solves for the bus voltage angles at all but the reference bus,
/// given the full system B matrix and the vector of bus real power injections,
/// the initial vector of bus voltage angles (in radians), and column vectors
/// with the lists of bus indices for the swing bus, PV buses, and PQ buses,
/// respectively. Returns a vector of bus voltage angles in radians.
pub(crate) fn dc_pf(
    b_mat: &CSR<usize, f64>,
    p_bus: &[f64],
    va0: &[f64],
    ref_: &[usize],
    pv: &[usize],
    pq: &[usize],
    solver: &dyn Solver<usize, f64>,
) -> Result<(Vec<f64>, bool)> {
    let va_threshold = 1e5; // arbitrary threshold on |Va| for declaring failure

    // initialize result vector
    let mut va = va0.to_vec();
    let mut success = true; // successful by default

    // update angles for non-reference buses
    let pvpq = [pv, pq].concat();

    // Va([pv; pq]) = B([pv; pq], [pv; pq]) \ ...
    //                     (Pbus([pv; pq]) - B([pv; pq], ref) * Va0(ref));

    let b_pvpq = b_mat.select(Some(&pvpq), Some(&pvpq))?;
    let b_ref = b_mat.select(Some(&pvpq), Some(ref_))?;
    let p_bus_pvpq: Vec<f64> = pvpq.iter().map(|&i| p_bus[i]).collect();
    let va_ref: Vec<f64> = ref_.iter().map(|&i| va0[i]).collect();

    let mut rhs: Vec<f64> = zip(p_bus_pvpq, b_ref * &va_ref)
        .map(|(p_bus, p_ref)| p_bus - p_ref)
        .collect();

    let va_pvpq = {
        solver.solve(
            b_pvpq.cols(),
            b_pvpq.colidx(),
            b_pvpq.rowptr(),
            b_pvpq.values(),
            &mut rhs,
            true,
        )?;
        rhs
    };

    pvpq.iter()
        .enumerate()
        .for_each(|(i, &j)| va[j] = va_pvpq[i]);
    // set_slice(&mut va, &pvpq, &va_pvpq);

    // if va.abs().max() > va_threshold {
    if norm_inf(&va).unwrap() > va_threshold {
        success = false;
    }

    Ok((va, success))
}
