//! This module takes care of event generation and storage

use ::{
    linalg::{
        Momentum,
        E,
        Matrix2x4,
        Matrix3,
        Matrix3x2,
        Matrix4x3,
        Matrix5,
        Matrix5x4,
        MatrixSlice,
        MatrixSliceMut,
        U1,
        U3,
        U4,
        U5,
        Vector2,
        Vector3,
        Vector5,
        X,
        xyz,
        Y,
        Z
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


/// Number of incoming particles
pub const INCOMING_COUNT: usize = 2;

/// Number of outgoing particles (replaces original INP)
pub const OUTGOING_COUNT: usize = 3;

/// Vectors and matrices sized to the amount of particles in an event (=5)
pub type ParticleVector<T> = Vector5<T>;
pub type ParticleMatrix<T> = Matrix5<T>;

/// Event data matrix definitions (columns are coordinates, rows are momenta)
type EventMatrix = Matrix5x4<Real>;
type OutgoingMomentaSlice<'a> = MatrixSlice<'a, Real, U3, U4, U1, U5>;
type OutgoingMomentaSliceMut<'a> = MatrixSliceMut<'a, Real, U3, U4, U1, U5>;

/// Row of the incoming electron in the event data matrix
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
    incoming_momenta: Matrix2x4<Real>,
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
        let incoming_momenta = Matrix2x4::new(
            -half_e_tot, 0., 0., half_e_tot,
            half_e_tot, 0., 0., half_e_tot
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
        // Pregenerate the random parameters to shield later computations from
        // the averse impact of RNG calls on the compiler's loop optimizations
        //
        // FIXME: This temporarily uses a different RNG order than 3photons.
        //        Should ultimately extract the fluctuation in a method.
        //
        assert_eq!(OUTGOING_COUNT, 3, "This code assumes 3 outgoing particles");
        let params = Matrix3::from_column_slice(&rng.random9()[..]);
        let cos_theta = params.fixed_columns::<U1>(0).map(|r| 2. * r - 1.);
        let exp_min_e = params.fixed_columns::<U1>(1)
                              .component_mul(&params.fixed_columns::<U1>(2));
        let sincos_phi_mat = Self::random_unit_3x2d(rng);

        // Generate massless outgoing 4-momenta in infinite phase space
        //
        // FIXME: Sincos coordinate order is reversed w.r.t. 3photons
        //
        // FIXME: The main obvious remaining bottleneck of this function is that
        //        it spends 25% of its time computing scalar logarithms. Using
        //        a vectorized ln() implementation should help there.
        //
        let sin_theta = cos_theta.map(|cos| sqrt(1. - sqr(cos)));
        let energy = exp_min_e.map(|e_me| -ln(e_me));
        let q_mat = Matrix4x3::from_fn(|coord, par| {
            energy[par] * match coord {
                X => sin_theta[par] * sincos_phi_mat[(par, X)],
                Y => sin_theta[par] * sincos_phi_mat[(par, Y)],
                Z => cos_theta[par],
                E => 1.,
                _ => unreachable!()
            }
        });

        // Calculate the parameters of the conformal transformation
        let r = &Momentum::from_fn(|coord, _| {
            q_mat.fixed_rows::<U1>(coord).iter().sum()
        });
        let r_norm_2 = r[E] * r[E] - xyz(r).norm_squared();
        let alpha = self.e_tot / r_norm_2;
        let r_norm = sqrt(r_norm_2);
        let beta = 1. / (r_norm + r[E]);

        // Perform the conformal transformation from Q's to output 4-momenta
        let tr_q_mat = q_mat.transpose();
        let rq = tr_q_mat.fixed_columns::<U3>(X) * xyz(r);
        let q_e = tr_q_mat.fixed_columns::<U1>(E);
        let p_e = alpha * (r[E] * q_e - rq);
        let tr_q_xyz = tr_q_mat.fixed_columns::<U3>(X);
        let b_rq_e = beta * rq - q_e;
        let p_xyz = alpha * (r_norm * tr_q_xyz + b_rq_e * xyz(r).transpose());

        // Build the final event: incoming momenta + output 4-momenta
        let mut event = Event(Matrix5x4::from_fn(|par, coord| {
            if par < INCOMING_COUNT {
                self.incoming_momenta[(par, coord)]
            } else if coord <= Z {
                p_xyz[(par-INCOMING_COUNT, coord)]
            } else if coord == E {
                p_e[par-INCOMING_COUNT]
            } else {
                unreachable!()
            }
        }));

        // Sort the output 4-momenta in order of decreasing energy
        if cfg!(not(feature = "no-photon-sorting")) {
            let mut outgoing = event.outgoing_momenta_mut();
            for par1 in 0..OUTGOING_COUNT-1 {
                for par2 in par1+1..OUTGOING_COUNT {
                    if outgoing[(par2, E)] > outgoing[(par1, E)] {
                        outgoing.swap_rows(par1, par2);
                    }
                }
            }
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
    /// FIXME: Uses a different RNG order w.r.t. original 3photons
    ///
    fn random_unit_3x2d(rng: &mut RandomGenerator) -> Matrix3x2<Real> {
        assert_eq!(OUTGOING_COUNT, 3, "This code assumes 3 outgoing particles");

        // This function has two operating modes: a default mode which produces
        // bitwise identical results w.r.t. the original 3photons code, and a
        // mode which uses a different (faster) algorithm.
        if cfg!(feature = "fast-sincos") {
            // Grab random points on the unit square
            let mut p_mat = Matrix3x2::from_iterator(
                rng.random6().iter().map(|r| 2. * r - 1.)
            );

            // Re-roll each point until it falls on the unit disc, and is not
            // too close to the origin (the latter improves numerical stability)
            let mut n2_vec = Vector3::from_fn(|part, _|
                p_mat.fixed_rows::<U1>(part).norm_squared()
            );
            for part in 0..OUTGOING_COUNT {
                const MIN_POSITIVE_2: Real = MIN_POSITIVE * MIN_POSITIVE;
                while n2_vec[part] > 1. || n2_vec[part] < MIN_POSITIVE_2 {
                    let new_row = Vector2::from_iterator(
                        rng.random2().iter().map(|r| 2. * r - 1.)
                    );
                    p_mat.set_row(part, &new_row.transpose());
                    n2_vec[part] = new_row.norm_squared();
                }
            }

            // Now you only need to normalize to get points on the unit circle
            let norm = n2_vec.map(|n2| 1. / sqrt(n2));
            for part in 0..OUTGOING_COUNT {
                p_mat.fixed_rows_mut::<U1>(part).apply(|c| c * norm[part]);
            }
            p_mat
        } else {
            // This code path strictly follows the original 3photons algorithm
            let phi = Vector3::from_iterator(
                rng.random3().iter().map(|r| 2. * PI * r)
            );
            Matrix3x2::from_fn(|par, coord| {
                match coord {
                    X => cos(phi[par]),
                    Y => sin(phi[par]),
                    _ => unreachable!()
                }
            })
        }
    }

    /// Simulate the impact of a certain number of calls to "generate()" on a
    /// random number generator. This code must be manually synchronized with
    /// the implementation of "generate()", but such is the price for perfect
    /// reproducibility between single-threaded and multi-threaded mode...
    ///
    /// FIXME: Update this if we stabilize the new RNG order
    ///
    #[cfg(all(feature = "multi-threading",
              not(feature = "faster-threading")))]
    pub(crate) fn simulate_event_batch(rng: &mut RandomGenerator,
                                       num_events: usize) {
        if cfg!(feature = "fast-sincos") {
            // If fast-sincos is enabled, the number of RNG calls per event is
            // nondeterministic, so we must simulate events one by one.
            for _ in 0..num_events*OUTGOING_COUNT {
                rng.skip(1);
                Self::random_unit_2d(rng);
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
///
/// Encapsulates a vector of incoming and outgoing 4-momenta
///
pub struct Event(EventMatrix);
//
impl Event {
    // ### ACCESSORS ###

    /// Access the full internal 4-momentum array by reference
    pub fn all_momenta(&self) -> &EventMatrix {
        &self.0
    }

    /// Access the momentum of a single particle (internal for now)
    fn momentum(&self, par: usize) -> Momentum {
        Momentum::from_iterator(self.0.fixed_rows::<U1>(par).iter().cloned())
    }

    /// Access the electron 4-momentum only
    pub fn electron_momentum(&self) -> Momentum {
        self.momentum(INCOMING_E_M)
    }

    /// Access the positron 4-momentum only
    #[allow(dead_code)]
    pub fn positron_momentum(&self) -> Momentum {
        self.momentum(INCOMING_E_P)
    }

    /// Access the momentum of a single outgoing photon
    pub fn outgoing_momentum(&self, par: usize) -> Momentum {
        self.momentum(OUTGOING_SHIFT + par)
    }

    /// Access the outgoing 4-momenta only
    pub fn outgoing_momenta(&self) -> OutgoingMomentaSlice {
        self.0.fixed_rows::<U3>(OUTGOING_SHIFT)
    }

    /// Mutable access to the outgoing 4-momenta (for internal use)
    fn outgoing_momenta_mut(&mut self) -> OutgoingMomentaSliceMut {
        self.0.fixed_rows_mut::<U3>(OUTGOING_SHIFT)
    }

    /// Minimal outgoing photon energy
    pub fn min_photon_energy(&self) -> Real {
        if cfg!(feature = "no-photon-sorting") {
            let first_out_e = self.outgoing_momentum(0)[E];
            self.outgoing_momenta()
                .fixed_columns::<U1>(E)
                .iter()
                .skip(1)
                .fold(first_out_e, |e1, &e2| if e1 < e2 { e1 } else { e2 })
        } else {
            // Use the fact that photons are sorted by decreasing energy
            self.outgoing_momentum(OUTGOING_COUNT-1)[E]
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