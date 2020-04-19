//! Physical couplings used for result computations

use crate::{
    config::Configuration,
    numeric::{floats::consts::PI, Float},
};

use prefix_num_ops::real::*;

/// Set of physical couplings
pub struct Couplings {
    /// Standard Model contribution electromagnetic coupling √(4𝜋𝛼)³
    pub g_a: Float,

    /// 𝛽₊ anomalous contribution electroweak coupling
    pub g_beta_p: Float,

    /// 𝛽₋ anomalous contribution electroweak coupling
    pub g_beta_m: Float,
}
//
impl Couplings {
    /// Fill in the parameters using data from the configuration file
    pub fn new(cfg: &Configuration) -> Self {
        let e2 = 4. * PI * cfg.alpha;
        let e2_z = 4. * PI * cfg.alpha_z;
        let cos2_weinberg = 1. - cfg.sin2_weinberg;
        let g_beta = -sqrt(e2_z / (4. * cos2_weinberg * cfg.sin2_weinberg)) / cfg.m_z0.powi(4);
        Couplings {
            g_a: -sqrt(e2).powi(3),
            g_beta_p: g_beta,
            g_beta_m: g_beta,
        }
    }
}
