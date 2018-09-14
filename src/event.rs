//! This module takes care of event generation and storage

use ::{
    linalg::{Momentum, E, X, xyz, xyz_mut, Y, Z},
    numeric::{
        functions::{cos, exp, ln, sin, sqr, sqrt},
        Real,
        reals::{
            MIN_POSITIVE,
            consts::{FRAC_PI_2, PI},
        },
    },
    random::RandomGenerator,
};

use num_traits::identities::Zero;

use std::cmp::Ordering;


/// Number of incoming particles
pub const INCOMING_COUNT: usize = 2;

/// Number of outgoing particles (replaces original INP)
pub const OUTGOING_COUNT: usize = 3;

/// Total number of particles in an event
pub const PARTICLE_COUNT: usize = INCOMING_COUNT+OUTGOING_COUNT;

/// Index of the incoming electron in the 4-momentum array
pub const INCOMING_E_M: usize = 0;

/// Index of the incoming positron in the 4-momentum array
pub const INCOMING_E_P: usize = 1;

/// Index of the first outgoing photon in the 4-momentum array
pub const OUTGOING_SHIFT: usize = INCOMING_COUNT;


/// Generator of ee -> ppp events
pub struct EventGenerator {
    /// Total center-of-mass energy of the collision
    e_tot: Real,

    /// Weight of generated events
    ev_weight: Real,
}
//
impl EventGenerator {
    // ### CONSTRUCTION ###

    /// Initialize event generation for a center-of-mass energy of e_tot.
    ///
    /// Combines former functionality of ppp constructor and IBEGIN-based lazy
    /// initialization from the original C++ 3photons code.
    ///
    pub fn new(e_tot: Real) -> Self {
        // Check on the number of particles. The check for N<101 is gone since
        // unlike the original RAMBO, we don't use arrays of hardcoded size.
        assert!(OUTGOING_COUNT > 1);

        // As currently written, this code only works for two incoming particles
        assert_eq!(INCOMING_COUNT, 2);

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

        // All generated events will have the same weight: pre-compute it
        let ln_weight = (2. * (OUTGOING_COUNT as Real) - 4.) * ln(e_tot) + z;
        assert!((ln_weight >= -180.) && (ln_weight <= 174.));
        let ev_weight = exp(ln_weight);
        
        // Construct and return the output data structure
        EventGenerator {
            e_tot,
            ev_weight,
        }
    }


    // ### EVENT GENERATION ###

    /// Use a highly specialized version of the RAMBO (RAndom Momenta
    /// Beautifully Organized) algorithm from S.D. Ellis, R. Kleiss and W.J.
    /// Stirling to generate the 4-momenta of the three outgoing photons.
    ///
    /// All events have the same weight, it can be queried via event_weight().
    ///
    /// The 4-momenta of output photons are sorted by decreasing energy.
    ///
    pub fn generate(&self, rng: &mut RandomGenerator) -> Event {
        // Prepare storage for the final event
        let half_e_tot = self.e_tot / 2.;
        let mut event = Event { p: [Momentum::zero(); PARTICLE_COUNT] };
        event.p[INCOMING_E_M] = Momentum::new(-half_e_tot, 0., 0., half_e_tot);
        event.p[INCOMING_E_P] = Momentum::new(half_e_tot, 0., 0., half_e_tot);

        // TODO: There might be room for more vectorization or layout optims.
        //       Here is where I ended up the last time I looked at it:
        //
        //       - There is an obvious vectorization opportunity in energy
        //         computations (see below).
        //       - The normalization in the middle is a sequential bottleneck.
        //       - A vectorized RNG would be useful here and there.
        //       - It's not that obvious what the layout of q_arr shoud be. I
        //         think the computation of r works much better with q_arr in
        //         its current layout, but the computation of p at the end could
        //         work better with q_arr in a transposed layout.
        //       - The ideal layout for p is mostly determined from clients, but
        //         there is also an incentive to compute p in the same
        //         orientation as q.

        // Pregenerate the random parameters to shield later computations from
        // the averse impact of RNG calls on the compiler's loop optimizations
        #[derive(Clone, Copy, Default)]
        struct RandomParameters {
            cos_theta: Real,
            sincos_phi: [Real; 2],
            exp_minus_e: Real,
        }
        let mut rand_params_arr = [RandomParameters::default(); OUTGOING_COUNT];
        for rand_params in rand_params_arr.iter_mut() {
            *rand_params = RandomParameters {
                cos_theta: 2. * rng.random() - 1.,
                sincos_phi: Self::random_sincos(rng),
                exp_minus_e: rng.random() * rng.random(),
            };
        }

        // Generate massless outgoing 4-momenta in infinite phase space
        //
        // FIXME: The main obvious remaining bottleneck of this function is that
        //        it spends 25% of its time computing scalar logarithms. Using
        //        a vectorized ln() implementation should help there.
        //
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
        let r: &Momentum = &q_arr.iter().sum();
        let r_norm_2 = r[E] * r[E] - xyz(r).norm_squared();
        let alpha = self.e_tot / r_norm_2;
        let r_norm = sqrt(r_norm_2);
        let beta = 1. / (r_norm + r[E]);

        // Transform the Q's conformally into output 4-momenta
        for (p, q) in event.outgoing_momenta_mut().iter_mut()
                                                  .zip(q_arr.iter())
        {
            let rq = xyz(r).dot(&xyz(q));
            let p_xyz = r_norm * xyz(q) + (beta * rq - q[E]) * xyz(r);
            xyz_mut(p).copy_from(&p_xyz);
            p[E] = r[E] * q[E] - rq;
            *p *= alpha;
        }

        // Sort the output 4-momenta in order of decreasing energy
        event.outgoing_momenta_mut()
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

        // Hand off the generated event
        event
    }

    /// Generate a (sin(x), cos(x)) pair where x is uniform in [0, 2*PI[
    ///
    /// TODO: Make this a full-blown unit vector generator and explore if we can
    ///       do better on the theta front.
    ///
    /// FIXME: I would like to use a vector type instead of an array here, but
    ///        nalgebra's operations cause an unacceptable performance hit.
    ///
    fn random_sincos(rng: &mut RandomGenerator) -> [Real; 2] {
        // This function has two operating modes: a default mode which produces
        // bitwise identical results w.r.t. the original 3photons code, and a
        // mode which uses a different (faster) algorithm.
        if cfg!(feature = "fast-sincos") {
            // In a nutshell, this path is faster because it favors cheap RNG
            // calls over expensive trigonometric functions
            const MIN_POSITIVE_2: Real = MIN_POSITIVE * MIN_POSITIVE;
            loop {
                // Grab a point on the unit square
                let x = 2. * rng.random() - 1.;
                let y = 2. * rng.random() - 1.;

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
            // This code path strictly follows the original 3photons algorithm
            let phi = 2. * PI * rng.random();
            [cos(phi), sin(phi)]
        }
    }

    /// Simulate the impact of a certain number of calls to "generate()" on a
    /// random number generator. This code must be manually synchronized with
    /// the implementation of "generate()", but such is the price for perfect
    /// reproducibility between single-threaded and multi-threaded mode...
    #[cfg(all(feature = "multi-threading",
              not(feature = "faster-threading")))]
    pub(crate) fn simulate_event_batch(rng: &mut RandomGenerator,
                                       num_events: usize) {
        if cfg!(feature = "fast-sincos") {
            // If fast-sincos is enabled, the number of RNG calls per event is
            // nondeterministic, so we must simulate events one by one.
            for _ in 0..num_events*OUTGOING_COUNT {
                rng.skip(1);
                Self::random_sincos(rng);
                rng.skip(2);
            }
        } else {
            // If fast-sincos is not enabled, we know exactly how many RNG calls
            // will be made per event, and we can let the RNG skip through the
            // events as quickly as it can.
            rng.skip(num_events*OUTGOING_COUNT*4);
        }
    }


    // ### EVENT PROPERTIES ###

    /// Access the event weight (identical for all generated events)
    pub fn event_weight(&self) -> Real {
        self.ev_weight
    }
}


/// Storage for ee -> ppp event data
pub struct Event {
    /// Array of incoming and outgoing 4-momenta
    ///
    /// TODO: Should use a matrix here (and review row/column layout)
    ///
    p: [Momentum; PARTICLE_COUNT],
}
//
impl Event {
    // ### ACCESSORS ###

    /// Access the full internal 4-momentum array by reference
    pub fn all_momenta(&self) -> &[Momentum; PARTICLE_COUNT] {
        &self.p
    }

    /// Access the electron 4-momentum only
    pub fn electron_momentum(&self) -> &Momentum {
        &self.p[INCOMING_E_M]
    }

    /// Access the positron 4-momentum only
    #[allow(dead_code)]
    pub fn positron_momentum(&self) -> &Momentum {
        &self.p[INCOMING_E_P]
    }

    /// Access the outgoing 4-momenta only
    pub fn outgoing_momenta(&self) -> &[Momentum; OUTGOING_COUNT] {
        array_ref![self.p, OUTGOING_SHIFT, OUTGOING_COUNT]
    }

    /// Mutable access to the outgoing 4-momenta (for internal use)
    fn outgoing_momenta_mut(&mut self) -> &mut [Momentum; OUTGOING_COUNT] {
        array_mut_ref![self.p, OUTGOING_SHIFT, OUTGOING_COUNT]
    }

    /// Minimal outgoing photon energy
    pub fn min_photon_energy(&self) -> Real {
        // Use the fact that photons are sorted in order of decreasing energy
        self.p[OUTGOING_SHIFT+OUTGOING_COUNT-1][E]
    }


    // ### DEBUGGING TOOLS ###

    /// Dump 4-momenta of the 3 outgoing photons
    #[allow(dead_code)]
    pub fn display(&self) {
        let p_out = self.outgoing_momenta();
        for ip in p_out.iter().enumerate() {
            let (i, p) = ip;
            println!("p{}: {}", i+1, p);
        }
        println!();
    }
}