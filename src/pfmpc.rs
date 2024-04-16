use powers::MPC;
use std::ops::{Deref, DerefMut};

#[derive(Clone, Default)]
pub struct PFMPC {
    mpc: MPC,
    pub(crate) success: Option<bool>,
    pub(crate) iterations: Option<usize>,
}

impl PFMPC {
    pub(crate) fn new(mpc: MPC) -> Self {
        Self {
            mpc,
            success: None,
            iterations: None,
        }
    }

    pub fn success(&self) -> bool {
        self.success.unwrap_or_default()
    }

    pub fn iterations(&self) -> usize {
        self.iterations.unwrap_or_default()
    }
}

impl Deref for PFMPC {
    type Target = MPC;

    fn deref(&self) -> &Self::Target {
        &self.mpc
    }
}

impl DerefMut for PFMPC {
    fn deref_mut(&mut self) -> &mut Self::Target {
        &mut self.mpc
    }
}
