use crate::pfopt::MPOpt;
use anyhow::Result;
use caseformat::{Branch, Bus, Gen};
use powers::MPC;

pub(crate) fn radial_pf(
    _base_mva: f64,
    _bus0: &[Bus],
    _gen0: &[Gen],
    _branch0: &[Branch],
    _mpopt: &MPOpt,
) -> Result<(MPC, bool, usize)> {
    unimplemented!("radial power flow")
}
