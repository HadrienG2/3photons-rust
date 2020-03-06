///! This module allows integrating simulation results across generated events
use crate::{
    config::Configuration,
    event::{NUM_OUTGOING, NUM_SPINS},
    linalg::{dimension::*, vecmat::*},
    matelems::{MEsContributions, MEsVector, A, B_M, B_P, NUM_MAT_ELEMS, R_MX},
    numeric::{functions::*, reals::consts::PI, Float},
    resfin::{FinalResults, PerSpinMEs},
};

use num_traits::Zero;

/// This struct will accumulate intermediary results during integration, and
/// ultimately compute the final results (see FinalResults below).
pub struct ResultsAccumulator<'cfg> {
    // ### RESULT ACCUMULATORS ###
    /// Number of integrated events
    selected_events: usize,

    /// Accumulated cross-section for each contribution
    spm2: MEsVector,

    /// Accumulated variance for each contribution
    vars: MEsVector,

    /// Impact of each contribution on the cross-section
    sigma_contribs: MEsVector,

    /// Accumulated total cross-section
    sigma: Float,

    /// Accumulated total variance
    variance: Float,

    // ### PHYSICAL CONSTANTS (CACHED FOR FINALIZATION) ###
    /// Configuration of the simulation
    cfg: &'cfg Configuration,

    /// Common factor, non-averaged over spins=1
    ///                  /(symmetry factor)
    ///                  *(conversion factor GeV^-2->pb)
    /// To average over spins, add :
    ///                  /(number of incoming helicities)
    fact_com: Float,

    /// Event weight, with total phase space normalization
    norm_weight: Float,

    /// Z° propagator
    propag: Float,

    /// ??? (Ask Vincent Lafage)
    ecart_pic: Float,
}
//
impl<'cfg> ResultsAccumulator<'cfg> {
    /// Prepare for results integration
    pub fn new(cfg: &'cfg Configuration, event_weight: Float) -> Self {
        // This code depends on some aspects of the problem definition
        assert_eq!(NUM_MAT_ELEMS, 5);

        // Common factor (see definition and remarks above)
        let fact_com = 1. / 6. * cfg.convers;
        let gzr = cfg.g_z0 / cfg.m_z0;

        // Sum over polarisations factors
        let p_aa = 2.;
        let p_ab = 1. - 4. * cfg.sin2_w;
        let p_bb = p_ab + 8. * sqr(cfg.sin2_w);

        // Homogeneity coefficient
        let c_aa = fact_com * p_aa;
        let c_ab = fact_com * p_ab / sqr(cfg.m_z0);
        let c_bb = fact_com * p_bb / powi(cfg.m_z0, 4);

        // Switch to dimensionless variable
        let dzeta = sqr(cfg.e_tot / cfg.m_z0);
        let ecart_pic = (dzeta - 1.) / gzr;
        let propag = 1. / (1. + sqr(ecart_pic));

        // Apply total phase space normalization to the event weight
        let n_ev = cfg.num_events as Float;
        let norm = powi(2. * PI, 4 - 3 * (NUM_OUTGOING as i32)) / n_ev;
        // NOTE: This replaces the original WTEV, previously reset every event
        let norm_weight = event_weight * norm;

        // Compute how much each result contribution adds to the cross-section.
        // Again, this avoids duplicate work in the integration loop.
        let com_contrib = norm_weight / 4.;
        let aa_contrib = com_contrib * c_aa;
        let bb_contrib = com_contrib * c_bb * propag / sqr(gzr);
        let ab_contrib = com_contrib * c_ab * 2. * cfg.beta_plus * propag / gzr;
        let sigma_contribs = MEsVector::new(
            aa_contrib,                       // A
            bb_contrib * sqr(cfg.beta_plus),  // B_P
            bb_contrib * sqr(cfg.beta_minus), // B_M
            ab_contrib * ecart_pic,           // R_MX
            -ab_contrib,                      // I_MX
        );

        // Return a complete results builder
        ResultsAccumulator {
            selected_events: 0,
            spm2: MEsVector::zero(),
            vars: MEsVector::zero(),
            sigma_contribs,
            sigma: 0.,
            variance: 0.,

            cfg,
            fact_com,
            norm_weight,
            propag,
            ecart_pic,
        }
    }

    /// Integrate one intermediary result into the simulation results
    #[allow(clippy::needless_pass_by_value)]
    pub fn integrate(&mut self, result: MEsContributions) {
        self.selected_events += 1;
        let spm2_dif = result.m2_sums();
        self.spm2 += spm2_dif;
        self.vars += spm2_dif.map(sqr);
        let weight = spm2_dif.dot(&self.sigma_contribs);
        self.sigma += weight;
        self.variance += sqr(weight);
    }

    /// Integrate simulation results from another ResultsAccumulator
    #[allow(clippy::needless_pass_by_value)]
    pub fn merge(&mut self, other: Self) {
        self.selected_events += other.selected_events;
        self.spm2 += other.spm2;
        self.vars += other.vars;
        self.sigma += other.sigma;
        self.variance += other.variance;
    }

    /// Turn integrated simulation data into finalized results
    pub fn finalize(mut self) -> FinalResults<'cfg> {
        // This code depends on some aspects of the problem definition
        assert_eq!(NUM_SPINS, 2);
        assert_eq!(NUM_MAT_ELEMS, 5);

        // Simulation configuration shorthand
        let cfg = self.cfg;

        // Keep around a floating-point version of the total event count
        let n_ev = cfg.num_events as Float;

        // Compute the relative uncertainties for one spin
        for (&v_spm2, v_var) in self.spm2.iter().zip(self.vars.iter_mut()) {
            *v_var = (*v_var - sqr(v_spm2) / n_ev) / (n_ev - 1.);
            *v_var = sqrt(*v_var / n_ev) / abs(v_spm2 / n_ev);
        }

        // Copy for the opposite spin
        let mut spm2 = PerSpinMEs::from_fn(|_spin, res| self.spm2[res]);
        let vars = PerSpinMEs::from_fn(|_spin, res| self.vars[res]);

        // Electroweak polarisations factors for the 𝛽₊/𝛽₋ anomalous
        // contribution
        let pol_p = -2. * cfg.sin2_w;
        let pol_m = 1. + pol_p;
        let pols = Vector2::new(pol_m, pol_p);

        // Take polarisations into account
        spm2.fixed_columns_mut::<U4>(B_P)
            .column_iter_mut()
            .for_each(|mut col| col.component_mul_assign(&pols));
        spm2.fixed_columns_mut::<U2>(B_P)
            .column_iter_mut()
            .for_each(|mut col| col.component_mul_assign(&pols));

        // Flux factor (=1/2s for 2 initial massless particles)
        let flux = 1. / (2. * sqr(cfg.e_tot));

        // Apply physical coefficients and Z⁰ propagator to each spin
        spm2 *= self.fact_com * flux * self.norm_weight;
        let gm_z0 = cfg.g_z0 * cfg.m_z0;
        spm2.fixed_columns_mut::<U4>(B_P)
            .apply(|x| x * self.propag / gm_z0);
        spm2.fixed_columns_mut::<U2>(B_P).apply(|x| x / gm_z0);
        spm2.column_mut(R_MX).apply(|x| x * self.ecart_pic);

        // Compute other parts of the result
        let beta_min = sqrt(spm2.column(A).sum() / spm2.column(B_P).sum());

        let ss_denom = spm2.column(A).sum();
        let ss_norm = 1. / (2. * sqrt(ss_denom));

        let ss_p = spm2.column(B_P).sum() * ss_norm;
        let ss_m = spm2.column(B_M).sum() * ss_norm;

        let inc_num = |col| spm2.column(col).component_mul(&vars.column(col)).norm();
        let inc_ss_common = inc_num(A) / (2. * abs(ss_denom));
        let inc = |col| inc_num(col) / abs(spm2.column(col).sum()) + inc_ss_common;
        let inc_ss_p = inc(B_P);
        let inc_ss_m = inc(B_M);

        let variance = (self.variance - sqr(self.sigma) / n_ev) / (n_ev - 1.);
        let prec = sqrt(variance / n_ev) / abs(self.sigma / n_ev);
        let sigma = self.sigma * flux;

        // Return the final results
        FinalResults {
            selected_events: self.selected_events,
            spm2,
            vars,
            sigma,
            variance,
            beta_min,
            prec,
            ss_p,
            inc_ss_p,
            ss_m,
            inc_ss_m,
            cfg,
        }
    }
}
