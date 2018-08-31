//! This module takes care of event generation and storage

use ::{
    linalg::{self, Momentum, X, Y, Z, E},
    numeric::{
        cos,
        exp,
        ln,
        Real,
        reals::{
            MIN_POSITIVE,
            consts::{FRAC_PI_2, PI},
        },
        sin,
        sqr,
        sqrt
    },
    ranf::RanfGenerator,
};

use num_traits::identities::Zero;

use std::cmp::Ordering;


/// Number of 4-impulsions that are generated per event (replaces original INP)
pub const OUTGOING_COUNT: usize = 3;

/// Incoming (electron/positron) and outgoing (photon) particle impulsions will
/// be stored together as consecutive rows in an array
pub const INCOMING_E_M: usize = 0;
pub const INCOMING_E_P: usize = 1;
pub const INCOMING_COUNT: usize = 2;
pub const OUTGOING_SHIFT: usize = INCOMING_COUNT;
pub const PARTICLE_COUNT: usize = INCOMING_COUNT+OUTGOING_COUNT;
type Momenta = [Momentum; PARTICLE_COUNT];

/// A client can get access to the outgoing momenta alone if needed.
/// The resulting slice will always have OUTGOING_COUNT elements.
pub type OutgoingMomenta<'a> = &'a [Momentum];


/// Array of incoming and generated particle momenta
pub struct EventGenerator {
    /// Total center-of-mass energy of the collision
    e_tot: Real,

    /// Weight of generated events
    ev_weight: Real,

    /// Random number generator
    rng: RanfGenerator,
}
//
impl EventGenerator {
    // ### CONSTRUCTION ###

    /// Initialize incoming 4-momenta with a center-of-mass energy of e_tot.
    /// Combines former functionality of ppp constructor and IBEGIN-based lazy
    /// initialization from the original C++ 3photons code.
    pub fn new(e_tot: Real) -> Self {
        // Check on the number of particles. The check for N<101 is gone since
        // unlike the original RAMBO, we don't use arrays of hardcoded size.
        debug_assert!(OUTGOING_COUNT > 1);

        // As currently written, this code only works for two incoming particles
        debug_assert_eq!(INCOMING_COUNT, 2);

        // Compute some numerical constants. Replaces the lazy initialization
        // from the original RAMBO code with something less branchy.
        println!("IBegin");
        let po2log = ln(FRAC_PI_2);
        let mut z = po2log;  // Replaces Z[INP-1] in the original 3photons code
        for k in 2..OUTGOING_COUNT {
            z += po2log - 2. * ln((k-1) as Real);
        }
        let z = z - ln((OUTGOING_COUNT-1) as Real);

        // NOTE: The check on total energy is gone, because we only generate
        //       massless photons and so the total energy will always be enough.
        //       Counting of nonzero masses is also gone because it was unused.

        // All generated events will have the same weight, so we pre-compute it
        let ln_weight = (2. * (OUTGOING_COUNT as Real) - 4.) * ln(e_tot) + z;
        assert!((ln_weight >= -180.) && (ln_weight <= 174.));
        let ev_weight = exp(ln_weight);
        
        // Construct and return the output data structure
        EventGenerator {
            e_tot,
            ev_weight,
            rng: RanfGenerator::new(),
        }
    }


    // ### EVENT GENERATION ###

    /// Use a highly specialized version of the RAMBO (RAndom Momenta
    /// Beautifully Organized) algorithm from S.D. Ellis, R. Kleiss and W.J.
    /// Stirling in order to generate the momenta of the three outgoing photons.
    ///
    /// All events have the same weight, which can be read via event_weight().
    ///
    /// The momenta of output photons is sorted in order of decreasing energy.
    ///
    pub fn generate(&mut self) -> Event {
        // Prepare storage for the final event
        let half_e_tot = self.e_tot / 2.;
        let mut event = Event { p: [Momentum::zero(); PARTICLE_COUNT] };
        event.p[INCOMING_E_M] = Momentum::new(-half_e_tot, 0., 0., half_e_tot);
        event.p[INCOMING_E_P] = Momentum::new(half_e_tot, 0., 0., half_e_tot);

        // Pregenerate the random parameters (thusly isolating the harmful
        // branching and data dependency side-effects of the RNG)
        #[derive(Clone, Copy, Default)]
        struct RandomParameters {
            cos_theta: Real,
            sincos_phi: [Real; 2],
            exp_minus_e: Real,
        }
        let mut rand_params_arr = [RandomParameters::default(); OUTGOING_COUNT];
        for rand_params in rand_params_arr.iter_mut() {
            *rand_params = RandomParameters {
                cos_theta: 2. * self.rng.random() - 1.,
                sincos_phi: self.random_sincos(),
                exp_minus_e: self.rng.random() * self.rng.random(),
            };
        }

        // Generate massless outgoing momenta in infinite phase space
        // TODO: Once Rust supports it, initialize q_arr more directly
        let mut q_arr = [Momentum::zero(); OUTGOING_COUNT];
        for (q, params) in q_arr.iter_mut().zip(rand_params_arr.iter()) {
            let cos_theta = params.cos_theta;
            let sin_theta = sqrt(1.0 - sqr(cos_theta));
            let [cos_phi, sin_phi] = params.sincos_phi;
            q[E] = - ln(params.exp_minus_e);
            q[Z] = q[E] * cos_theta;
            q[Y] = q[E] * sin_theta * cos_phi;
            q[X] = q[E] * sin_theta * sin_phi;
        }
        let q_arr = q_arr;

        // Calculate the parameters of the conformal transformation
        let r = q_arr.iter().sum();
        let r_xyz = linalg::xyz(&r);
        let r_norm_2 = r[E] * r[E] - r_xyz.dot(&r_xyz);
        let alpha = self.e_tot / r_norm_2;
        let r_norm = sqrt(r_norm_2);
        let beta = 1. / (r_norm + r[E]);

        // Transform the Q's conformally into output momenta
        for (p, q) in event.dump_outgoing_mut().iter_mut()
                                               .zip(q_arr.iter())
        {
            let q_xyz = linalg::xyz(q);
            let rq = r_xyz.dot(&q_xyz);
            let p_xyz = r_norm * q_xyz + (beta * rq - q[E]) * r_xyz;
            linalg::xyz_mut(p).copy_from(&p_xyz);
            p[E] = r[E] * q[E] - rq;
            *p *= alpha;
        }

        // Sort the output momenta in order of decreasing energy
        event.sort_output_momenta();

        // Hand off the generated event
        event
    }

    /// Generate a (sin(x), cos(x)) pair where x is uniform in [0, 2*PI[
    fn random_sincos(&mut self) -> [Real; 2] {
        // This function has two operating modes: a default mode which produces
        // bitwise identical results w.r.t. the original 3photons code, and a
        // mode which uses a different algorithm to go faster.
        if cfg!(feature = "fast-sincos") {
            // This code path takes advantage of the fact that random
            // number generation is much faster than sin/cos calls in
            // order to speed things up.
            const MIN_POSITIVE_2: Real = MIN_POSITIVE * MIN_POSITIVE;
            loop {
                // Grab a point on the unit square
                let x = 2. * self.rng.random() - 1.;
                let y = 2. * self.rng.random() - 1.;

                // Compute (squared) distance from center
                let n2 = sqr(x) + sqr(y);

                // Discard points which are outside of the unit circle
                // or too close to the center for good normalization.
                if n2 <= 1. && n2 >= MIN_POSITIVE_2 {
                    // Normalize by n and you get a point on the unit
                    // circle, i.e. a sin/cos pair!
                    let n = sqrt(n2);
                    break [x/n, y/n];
                }
            }
        } else {
            // This code path strictly follows the original 3photons alg
            let phi = 2. * PI * self.rng.random();
            [cos(phi), sin(phi)]
        }
    }


    // ### EVENT PROPERTIES ###

    /// Access the event weight (identical for all generated events)
    pub fn event_weight(&self) -> Real {
        self.ev_weight
    }
}


/// Array of incoming and generated particle momenta
pub struct Event {
    /// Array of incoming and outgoing 4-momenta
    p: Momenta,
}
//
impl Event {
    // ### ACCESSORS ###

    /// Access the full internal momentum matrix by reference
    pub fn dump_momenta(&self) -> &Momenta {
        &self.p
    }

    /// Access the electron momentum only
    pub fn dump_electron(&self) -> &Momentum {
        &self.p[INCOMING_E_M]
    }

    /// Access the positron momentum only
    #[allow(dead_code)]
    pub fn dump_positron(&self) -> &Momentum {
        &self.p[INCOMING_E_P]
    }

    /// Access the outgoing momenta only
    pub fn dump_outgoing(&self) -> OutgoingMomenta {
        &self.p[OUTGOING_SHIFT..PARTICLE_COUNT]
    }

    /// Mutable access to the outgoing momenta (not publicly exported)
    fn dump_outgoing_mut(&mut self) -> &mut [Momentum] {
        &mut self.p[OUTGOING_SHIFT..PARTICLE_COUNT]
    }


    // ### GENERATOR HELPERS ###

    /// Sort outgoing photons in order of decreasing energy
    /// Roughly equivalent to the original ppp::TRI method
    fn sort_output_momenta(&mut self) {
        self.dump_outgoing_mut()
            .sort_unstable_by(|a, b| {
                // Treat NaNs as equal
                if a[E] > b[E] {
                    Ordering::Less
                } else if a[E] < b[E] {
                    Ordering::Greater
                } else {
                    Ordering::Equal
                }
            });
    }


    // ### DEBUGGING TOOLS ###

    /// Dump 4-momenta of the 3 outgoing photons
    #[allow(dead_code)]
    pub fn display(&self) {
        let p_out = self.dump_outgoing();
        for ip in p_out.iter().enumerate() {
            let (i, p) = ip;
            println!("p{}: {}", i+1, p);
        }
        println!();
    }
}