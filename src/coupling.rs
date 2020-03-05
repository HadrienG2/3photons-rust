//! Physical couplings used for result computations

use crate::{
    config::Configuration,
    numeric::{functions::*, reals::consts::PI, Float},
};

/// Set of physical couplings
pub struct Couplings {
    /// Standard Model contribution electromagnetic coupling √(4𝜋𝛼)³
    pub g_a: Float,

    /// 𝛽₊ anomalous contribution electroweak coupling
    pub g_bp: Float,

    /// 𝛽₋ anomalous contribution electroweak coupling
    pub g_bm: Float,
}
//
impl Couplings {
    /// Fill in the parameters using data from the configuration file
    pub fn new(cfg: &Configuration) -> Self {
        let e2 = 4. * PI * cfg.alpha;
        let e2_z = 4. * PI * cfg.alpha_z;
        let cos2_w = 1. - cfg.sin2_w;
        let g_beta = -sqrt(e2_z / (4. * cos2_w * cfg.sin2_w)) / powi(cfg.m_z0, 4);
        Couplings {
            g_a: -powi(sqrt(e2), 3),
            g_bp: g_beta,
            g_bm: g_beta,
        }
    }
}
