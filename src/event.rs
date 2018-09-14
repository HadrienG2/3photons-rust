//! This module takes care of event generation and storage

use ::{
    linalg::{
        Momentum,
        E,
        Matrix5,
        U3,
        U5,
        Vector2,
        Vector3,
        Vector5,
        VectorSlice,
        VectorSliceMut,
        X,
        xyz,
        xyz_mut,
        Y
    },
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


/// Number of incoming particles
pub const INCOMING_COUNT: usize = 2;
type IncomingVector<T> = Vector2<T>;

/// Number of outgoing particles (replaces original INP)
pub const OUTGOING_COUNT: usize = 3;
type OutgoingVector<T> = Vector3<T>;
pub type OutgoingVectorSlice<'a, T> = VectorSlice<'a, T, U3, U5>;
type OutgoingVectorSliceMut<'a, T> = VectorSliceMut<'a, T, U3, U5>;

/// Total number of particles in an event (= sum of the above)
pub type ParticleVector<T> = Vector5<T>;
pub type ParticleMatrix<T> = Matrix5<T>;

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

    /// Incoming electron and positron momenta
    ///
    /// FIXME: Should be a matrix
    ///
    incoming_momenta: IncomingVector<Momentum>,
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

        // Compute the incoming particle momenta
        let half_e_tot = e_tot / 2.;
        let incoming_momenta = IncomingVector::new(
            Momentum::new(-half_e_tot, 0., 0., half_e_tot),
            Momentum::new(half_e_tot, 0., 0., half_e_tot),
        );
        
        // Construct and return the output data structure
        EventGenerator {
            e_tot,
            ev_weight,
            incoming_momenta,
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
        let mut event = Event { p: ParticleVector::zero() };
        event.p[INCOMING_E_M] = self.incoming_momenta[INCOMING_E_M];
        event.p[INCOMING_E_P] = self.incoming_momenta[INCOMING_E_P];

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
        const COS_THETA: usize = 0;
        const COS_PHI: usize = 1;
        const SIN_PHI: usize = 2;
        const EXP_MINUS_E: usize = 3;
        let rand_params_arr = OutgoingVector::from_fn(|_, _| {
            let cos_theta = 2. * rng.random() - 1.;
            let sincos_phi = Self::random_unit_2d(rng);
            let exp_minus_e = rng.random() * rng.random();
            [cos_theta, sincos_phi[X], sincos_phi[Y], exp_minus_e]
        });

        // Generate massless outgoing 4-momenta in infinite phase space
        //
        // FIXME: The main obvious remaining bottleneck of this function is that
        //        it spends 25% of its time computing scalar logarithms. Using
        //        a vectorized ln() implementation should help there.
        //
        let q_arr = OutgoingVector::from_iterator(
            rand_params_arr.into_iter().map(|rand_params| {
                let cos_theta = rand_params[COS_THETA];
                let sin_theta = sqrt(1. - sqr(cos_theta));
                let energy = -ln(rand_params[EXP_MINUS_E]);
                energy * Momentum::new(sin_theta * rand_params[SIN_PHI],
                                       sin_theta * rand_params[COS_PHI],
                                       cos_theta,
                                       1.)
            })
        );

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
        if cfg!(not(feature = "no-photon-sorting")) {
            assert_eq!(OUTGOING_COUNT, 3,
                       "This code assumes that there are 3 outgoing particles");
            let mut outgoing = event.outgoing_momenta_mut();
            if outgoing[1][E] > outgoing[0][E] { outgoing.swap_rows(0, 1); }
            if outgoing[2][E] > outgoing[0][E] { outgoing.swap_rows(0, 2); }
            if outgoing[2][E] > outgoing[1][E] { outgoing.swap_rows(1, 2); }
        }

        // Hand off the generated event
        event
    }

    /// Generate a vector on the unit circle with uniform angle distribution
    ///
    /// NOTE: Similar techniques may be used to generate a vector on the unit
    ///       sphere, but that benchmarked unfavorably, likely because it
    ///       entails bringing more computations close to the RNG calls and
    ///       because the 2D case fits available vector hardware more tightly.
    ///
    fn random_unit_2d(rng: &mut RandomGenerator) -> Vector2<Real> {
        // This function has two operating modes: a default mode which produces
        // bitwise identical results w.r.t. the original 3photons code, and a
        // mode which uses a different (faster) algorithm.
        if cfg!(feature = "fast-sincos") {
            // In a nutshell, this path is faster because it favors cheap RNG
            // calls over expensive trigonometric functions
            const MIN_POSITIVE_2: Real = MIN_POSITIVE * MIN_POSITIVE;
            loop {
                // Grab a random point on the unit square
                let mut p = Vector2::from_fn(|_, _| 2. * rng.random() - 1.);

                // Compute (squared) distance from center
                let n2 = p.norm_squared();

                // Discard points which are outside of the unit circle
                // or too close to the center for good normalization.
                if n2 <= 1. && n2 >= MIN_POSITIVE_2 {
                    // Normalize by n and you get a point on the unit
                    // circle, i.e. a sin/cos pair!
                    p /= sqrt(n2);
                    break p;
                }
            }
        } else {
            // This code path strictly follows the original 3photons algorithm
            let phi = 2. * PI * rng.random();
            Vector2::new(cos(phi), sin(phi))
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
    /// TODO: Should use a matrix here (and review row/column layout), but won't
    ///       yet as nalgebra performance and ergonomics make it a bad trade-off
    ///
    p: ParticleVector<Momentum>,
}
//
impl Event {
    // ### ACCESSORS ###

    /// Access the full internal 4-momentum array by reference
    pub fn all_momenta(&self) -> &ParticleVector<Momentum> {
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
    pub fn outgoing_momenta(&self) -> OutgoingVectorSlice<Momentum> {
        self.p.fixed_rows::<U3>(OUTGOING_SHIFT)
    }

    /// Mutable access to the outgoing 4-momenta (for internal use)
    fn outgoing_momenta_mut(&mut self) -> OutgoingVectorSliceMut<Momentum> {
        self.p.fixed_rows_mut::<U3>(OUTGOING_SHIFT)
    }

    /// Minimal outgoing photon energy
    pub fn min_photon_energy(&self) -> Real {
        if cfg!(feature = "no-photon-sorting") {
            self.outgoing_momenta()
                .iter()
                .map(|p| p[E])
                .fold(self.p[0][E], |e1, e2| if e1 < e2 { e1 } else { e2 })
        } else {
            // Use the fact that photons are sorted by decreasing energy
            self.outgoing_momenta()[OUTGOING_COUNT-1][E]
        }
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