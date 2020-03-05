//! This module contains everything that is needed to compute, store, and
//! analyze the final results: differential cross-section, sum & variance

use crate::{
    config::Configuration,
    event::{NUM_OUTGOING, NUM_SPINS},
    linalg::{dimension::*, vecmat::*},
    matelems::{MEsContributions, MEsVector, A, B_M, B_P, I_MX, NUM_MAT_ELEMS, R_MX},
    numeric::{functions::*, reals::consts::PI, Complex, Float},
};

use num_traits::Zero;

/// Matrix of per-spin result contributions
///
/// Rows are spins, columns are result contributions (in the rescont.rs sense)
///
pub type PerSpinResults = Matrix2x5<Float>;

/// Index of negative spin data
pub const SP_M: usize = 0;

/// Index of positive spin data
pub const SP_P: usize = 1;

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
        let mut spm2 = PerSpinResults::from_fn(|_spin, res| self.spm2[res]);
        let vars = PerSpinResults::from_fn(|_spin, res| self.vars[res]);

        // Electroweak polarisations factors for the 𝛽₊/𝛽₋ anomalous
        // contribution
        let pol_p = -2. * cfg.sin2_w;
        let pol_m = 1. + pol_p;

        // Take polarisations into account
        // FIXME: Revisit those ugly slices once const generics allow for it
        spm2.fixed_slice_mut::<U1, U2>(SP_M, B_P)
            .apply(|x| x * sqr(pol_m));
        spm2.fixed_slice_mut::<U1, U2>(SP_P, B_P)
            .apply(|x| x * sqr(pol_p));
        spm2.fixed_slice_mut::<U1, U2>(SP_M, R_MX)
            .apply(|x| x * pol_m);
        spm2.fixed_slice_mut::<U1, U2>(SP_P, R_MX)
            .apply(|x| x * pol_p);

        // Flux factor (=1/2s for 2 initial massless particles)
        let flux = 1. / (2. * sqr(cfg.e_tot));

        // Apply physical coefficients and Z⁰ propagator to each spin
        spm2 *= self.fact_com * flux * self.norm_weight;
        let gm_z0 = cfg.g_z0 * cfg.m_z0;
        spm2.fixed_columns_mut::<U4>(B_P)
            .apply(|x| x * self.propag / gm_z0);
        spm2.fixed_columns_mut::<U2>(B_P).apply(|x| x / gm_z0);
        spm2.fixed_columns_mut::<U1>(R_MX)
            .apply(|x| x * self.ecart_pic);

        // Compute other parts of the result
        // FIXME: Vectorize over spins once const generics make it less ugly
        let beta_min =
            sqrt((spm2[(SP_M, A)] + spm2[(SP_P, A)]) / (spm2[(SP_M, B_P)] + spm2[(SP_P, B_P)]));

        let ss_denom = spm2[(SP_M, A)] + spm2[(SP_P, A)];
        let ss_norm = 1. / (2. * sqrt(ss_denom));

        let ss_p = (spm2[(SP_M, B_P)] + spm2[(SP_P, B_P)]) * ss_norm;
        let ss_m = (spm2[(SP_M, B_M)] + spm2[(SP_P, B_M)]) * ss_norm;

        let inc_ss_common =
            sqrt(sqr(spm2[(SP_M, A)] * vars[(SP_M, A)]) + sqr(spm2[(SP_P, A)] * vars[(SP_P, A)]))
                / (2. * abs(ss_denom));

        let inc_ss_p = sqrt(
            sqr(spm2[(SP_M, B_P)] * vars[(SP_M, B_P)]) + sqr(spm2[(SP_P, B_P)] * vars[(SP_P, B_P)]),
        ) / abs(spm2[(SP_M, B_P)] + spm2[(SP_P, B_P)])
            + inc_ss_common;
        let inc_ss_m = sqrt(
            sqr(spm2[(SP_M, B_M)] * vars[(SP_M, B_M)]) + sqr(spm2[(SP_P, B_M)] * vars[(SP_P, B_M)]),
        ) / abs(spm2[(SP_M, B_M)] + spm2[(SP_P, B_M)])
            + inc_ss_common;

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

/// Final results of the simulation
pub struct FinalResults<'cfg> {
    /// Number of integrated events
    pub selected_events: usize,

    /// Cross-section for each spin
    pub spm2: PerSpinResults,

    /// Variance for each spin
    pub vars: PerSpinResults,

    /// Total cross-section
    pub sigma: Float,

    /// Relative precision
    pub prec: Float,

    /// Total variance
    pub variance: Float,

    /// Beta minimum (???)
    pub beta_min: Float,

    /// Statistical significance B+(pb-1/2) (???)
    pub ss_p: Float,

    /// Incertitide associated with ss_p
    pub inc_ss_p: Float,

    /// Statistical significance B-(pb-1/2) (???)
    pub ss_m: Float,

    /// Incertitude associated with ss_m
    pub inc_ss_m: Float,

    /// Configuration of the simulation (for further derivation)
    cfg: &'cfg Configuration,
}
//
impl<'cfg> FinalResults<'cfg> {
    /// Display results using Eric's (???) parametrization
    pub fn eric(&self) {
        assert_eq!(NUM_SPINS, 2);
        assert_eq!(NUM_MAT_ELEMS, 5);

        let cfg = self.cfg;

        // FIXME: Vectorize over spins once const generics make it less ugly
        let mu_th = cfg.br_ep_em * cfg.convers / (8. * 9. * 5. * sqr(PI) * cfg.m_z0 * cfg.g_z0);
        let lambda0_m = (-self.spm2[(SP_M, B_P)] + self.spm2[(SP_M, B_M)]) / 2.;
        let lambda0_p = (-self.spm2[(SP_P, B_P)] + self.spm2[(SP_P, B_M)]) / 2.;
        let mu0_m = (self.spm2[(SP_M, B_P)] + self.spm2[(SP_M, B_M)]) / 2.;
        let mu0_p = (self.spm2[(SP_P, B_P)] + self.spm2[(SP_P, B_M)]) / 2.;
        let mu_num = (self.spm2[(SP_M, B_P)]
            + self.spm2[(SP_M, B_M)]
            + self.spm2[(SP_P, B_P)]
            + self.spm2[(SP_P, B_M)])
            / 4.;

        println!();
        println!("       :        -          +");
        println!(
            "sigma0  : {:.6} | {:.6}",
            self.spm2[(SP_M, A)] / 2.,
            self.spm2[(SP_P, A)] / 2.
        );
        println!(
            "alpha0  : {:.5e} | {:.4e}",
            self.spm2[(SP_M, I_MX)] / 2.,
            self.spm2[(SP_P, I_MX)] / 2.
        );
        println!(
            "beta0   : {:} | {:}",
            -self.spm2[(SP_M, R_MX)] / 2.,
            -self.spm2[(SP_P, R_MX)] / 2.
        );
        println!("lambda0 : {:.4} | {:.4}", lambda0_m, lambda0_p);
        println!("mu0     : {:.4} | {:.5}", mu0_m, mu0_p);
        println!(
            "mu/lamb : {:.5} | {:.5}",
            mu0_m / lambda0_m,
            mu0_p / lambda0_p
        );
        println!("mu (num): {:.4}", mu_num);
        println!("rapport : {:.6}", mu_num / mu_th);
        println!("mu (th) : {:.4}", mu_th);
    }

    /// Display Fawzi's (???) analytical results and compare them to the Monte
    /// Carlo results that we have computed
    pub fn fawzi(&self) {
        assert_eq!(NUM_SPINS, 2);
        assert_eq!(NUM_MAT_ELEMS, 5);

        let cfg = self.cfg;
        let ev_cut = &cfg.event_cut;

        let mre = cfg.m_z0 / cfg.e_tot;
        let gre = cfg.g_z0 * cfg.m_z0 / sqr(cfg.e_tot);
        let x = 1. - sqr(mre);
        let sdz = Complex::new(x, -gre) / (sqr(x) + sqr(gre));
        let del = (1. - ev_cut.b_cut) / 2.;
        let eps = 2. * ev_cut.e_min / cfg.e_tot;
        let bra = cfg.m_z0 / (3. * 6. * powi(PI, 3) * 16. * 120.);
        let sig = 12. * PI / sqr(cfg.m_z0) * cfg.br_ep_em * cfg.g_z0 * bra / sqr(cfg.e_tot)
            * powi(cfg.e_tot / cfg.m_z0, 8)
            * sdz.norm_sqr()
            * cfg.convers;

        let eps_4 = powi(eps, 4);
        let del_2 = powi(del, 2);
        let del_3 = powi(del, 3);
        let f1 = 1. - 15. * eps_4 - 9. / 7. * (1. - 70. * eps_4) * del_2
            + 6. / 7. * (1. + 70. * eps_4) * del_3;
        let g1 = 1.
            - 30. * eps_4
            - 9. / 7. * (1. - 70. * eps_4) * del
            - 90. * eps_4 * del_2
            - 1. / 7. * (1. - 420. * eps_4) * del_3;
        let g2 = 1.
            - 25. * eps_4
            - 6. / 7. * (1. - 70. * eps_4) * del
            - 3. / 7. * (1. + 210. * eps_4) * del_2
            - 8. / 21. * (1. - 52.5 * eps_4) * del_3;
        let g3 = 1.
            - 195. / 11. * eps_4
            - 18. / 77. * (1. - 7. * eps_4) * del
            - 9. / 11. * (9. / 7. - 70. * eps_4) * del_2
            - 8. / 11. * (1. - 105. / 11. * eps_4) * del_3;

        let sincut_3 = powi(ev_cut.sin_cut, 3);
        let ff = f1 * (1. - sincut_3);
        let gg = g1 - 27. / 16. * g2 * ev_cut.sin_cut + 11. / 16. * g3 * sincut_3;

        let sig_p = sig * (ff + 2. * gg);
        let sig_m = sig_p + 2. * sig * gg;

        // FIXME: Vectorize over spins once const generics make it less ugly
        let mc_p = (self.spm2[(SP_M, B_P)] + self.spm2[(SP_P, B_P)]) / 4.;
        let mc_m = (self.spm2[(SP_M, B_M)] + self.spm2[(SP_P, B_M)]) / 4.;
        let incr_p = sqrt(
            sqr(self.spm2[(SP_M, B_P)] * self.vars[(SP_M, B_P)])
                + sqr(self.spm2[(SP_P, B_P)] * self.vars[(SP_P, B_P)]),
        ) / abs(self.spm2[(SP_M, B_P)] + self.spm2[(SP_P, B_P)]);
        let incr_m = sqrt(
            sqr(self.spm2[(SP_M, B_M)] * self.vars[(SP_M, B_M)])
                + sqr(self.spm2[(SP_P, B_M)] * self.vars[(SP_P, B_M)]),
        ) / abs(self.spm2[(SP_M, B_M)] + self.spm2[(SP_P, B_M)]);

        println!();
        println!("s (pb) :   Sig_cut_Th    Sig_Th      Rapport");
        println!("       :   Sig_Num");
        println!("       :   Ecart_relatif  Incertitude");
        println!();
        println!(
            "s+(pb) : {:.5} | {:.5} | {:.6}",
            sig_p,
            sig * 3.,
            sig_p / (3. * sig)
        );
        println!("       : {:.5}", mc_p);
        println!(
            "       : {:.6} | {:.8} | {:.2}",
            mc_p / sig_p - 1.,
            incr_p,
            (mc_p / sig_p - 1.) / incr_p
        );
        println!();
        println!(
            "s-(pb) : {:.5} | {:.4} | {:.6}",
            sig_m,
            sig * 5.,
            sig_m / (5. * sig)
        );
        println!("       : {:.5}", mc_m);
        println!(
            "       : {:.6} | {:.9} | {:.2}",
            mc_m / sig_m - 1.,
            incr_m,
            (mc_m / sig_m - 1.) / incr_m
        );
        println!();
    }
}
