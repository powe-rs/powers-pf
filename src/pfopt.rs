use clap::ValueEnum;

#[derive(Debug, PartialEq, Copy, Clone, ValueEnum)]
pub enum Alg {
    /// Newton's method.
    NR = 0,
    /// Fast-Decoupled method (BX version).
    FDBX = 1,
    /// Fast-Decoupled method (XB version).
    FDXB = 2,
    /// Gauss-Seidel method.
    GS = 3,
    /// Power/current/admittance summation method (radial networks only).
    SUM = 4,
}

pub enum Sum {
    POWER,
    CURRENT,
    ADMITTANCE,
}

/// Type of nodal balance equation.
pub enum NodalBalance {
    POWER,
    CURRENT,
}

pub enum BusVoltage {
    /// bus voltage variables represented in polar coordinates
    POLAR,
    /// bus voltage variables represented in cartesian coordinates
    CARTESIAN,
    /// Polar updates computed via modified cartesian Jacobian
    HYBRID,
}

#[derive(PartialEq, Copy, Clone)]
pub enum GenQLimits {
    IgnoreLimits = 0,
    // Simultaneous bus type conversion.
    Simultaneous = 1,
    // One-at-a-time bus type conversion.
    OneAtATime = 2,
}

#[derive(Default)]
pub struct MPOpt {
    // Linearized DC power flow that assumes lossless branches,
    // 1pu voltages and small voltage angle differences.
    pub dc: bool,

    pub pf: PFOpt,
    pub exp: ExpOpt,
}

pub struct PFOpt {
    // AC power flow algorithm.
    pub algorithm: Alg,

    // Termination tolerance on per unit P & Q mismatch. Default value is 1e-8.
    pub tolerance: f64,

    // Maximum number of iterations for Newton's method. Default value is 10.
    pub max_it_nr: usize,
    // Maximum number of iterations for fast decoupled method. Default value is 30.
    pub max_it_fd: usize,
    // Maximum number of iterations for Gauss-Seidel method. Default value is 1000.
    pub max_it_gs: usize,

    // Enforce gen reactive power limits at expense of |V|.
    pub enforce_q_limits: GenQLimits,

    pub current_balance: NodalBalance,
    pub v_cartesian: BusVoltage,
    pub summation_method: Sum,
}

impl Default for PFOpt {
    fn default() -> Self {
        Self {
            algorithm: Alg::NR,
            tolerance: 1e-8,
            max_it_nr: 10,
            max_it_fd: 30,
            max_it_gs: 1000,
            enforce_q_limits: GenQLimits::IgnoreLimits,
            current_balance: NodalBalance::POWER,
            v_cartesian: BusVoltage::POLAR,
            summation_method: Sum::POWER,
        }
    }
}

#[derive(Default)]
pub struct ExpOpt {
    pub sys_wide_zip_loads: SysWideZipLoads,
}

#[derive(Default)]
pub struct SysWideZipLoads {
    pub pw: Option<[f64; 3]>,
    pub qw: Option<[f64; 3]>,
}
