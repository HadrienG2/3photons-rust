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
        U1,
        U3,
        U4,
        U5,
        Vector2,
        Vector3,
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

use num_traits::Zero;


/// Number of incoming particles
pub const INCOMING_COUNT: usize = 2;

/// Number of outgoing particles (replaces original INP)
pub const OUTGOING_COUNT: usize = 3;

/// Square matrix dimensioned to the amount of particles in an event (=5)
pub type ParticleMatrix<T> = Matrix5<T>;

/// Event data matrix definitions (columns are coordinates, rows are particles)
type EventMatrix = Matrix5x4<Real>;
type OutgoingMomentaSlice<'a> = MatrixSlice<'a, Real, U3, U4, U1, U5>;

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
        // Generate massless outgoing 4-momenta in infinite phase space
        let q = Self::generate_raw(rng);

        // Calculate the parameters of the conformal transformation
        let r = &Momentum::from_fn(|coord, _| {
            q.fixed_rows::<U1>(coord).iter().sum()
        });
        let r_norm_2 = r[E] * r[E] - xyz(r).norm_squared();
        let alpha = self.e_tot / r_norm_2;
        let r_norm = sqrt(r_norm_2);
        let beta = 1. / (r_norm + r[E]);

        // Perform the conformal transformation from Q's to output 4-momenta
        let tr_q = q.transpose();
        let rq = tr_q.fixed_columns::<U3>(X) * xyz(r);
        let q_e = tr_q.fixed_columns::<U1>(E);
        let mut p_e = alpha * (r[E]*q_e - rq);
        let tr_q_xyz = tr_q.fixed_columns::<U3>(X);
        let b_rq_e = beta * rq - q_e;
        let mut p_xyz = alpha * (r_norm*tr_q_xyz + b_rq_e*xyz(r).transpose());

        // Sort the output 4-momenta in order of decreasing energy (if enabled)
        if cfg!(not(feature = "no-photon-sorting")) {
            for par1 in 0..OUTGOING_COUNT-1 {
                for par2 in par1+1..OUTGOING_COUNT {
                    if p_e[par2] > p_e[par1] {
                        p_e.swap_rows(par1, par2);
                        // FIXME: swap_rows triggers UB here, figure out why
                        let mom1 = p_xyz.fixed_rows::<U1>(par1).into_owned();
                        let mom2 = p_xyz.fixed_rows::<U1>(par2).into_owned();
                        p_xyz.fixed_rows_mut::<U1>(par1).copy_from(&mom2);
                        p_xyz.fixed_rows_mut::<U1>(par2).copy_from(&mom1);
                    }
                }
            }
        }

        // Build the final event: incoming momenta + output 4-momenta
        let event = Event(Matrix5x4::from_fn(|par, coord| {
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

        // Hand off the generated event
        event
    }

    /// Generate massless outgoing 4-momenta in infinite phase space
    ///
    /// The output momenta are provided as a matrix where rows are momentum
    /// components and columns are particles.
    ///
    /// FIXME: The main obvious remaining bottleneck of this function is that
    ///        it spends 25% of its time computing scalar logarithms. Using
    ///        a vectorized ln() implementation should help there.
    ///
    fn generate_raw(rng: &mut RandomGenerator) -> Matrix4x3<Real> {
        assert_eq!(OUTGOING_COUNT, 3, "This code assumes 3 outgoing particles");

        // Pregenerate the random parameters to shield later computations from
        // the averse impact of RNG calls on the compiler's loop optimizations
        let mut cos_theta;
        let mut exp_min_e;
        let mut sincos_phi;
        if cfg!(feature = "simd-friendly-rng") {
            // Performance-oriented version
            let params = Matrix3::from_column_slice(&rng.random9()[..]);
            cos_theta = params.fixed_columns::<U1>(0).map(|r| 2. * r - 1.);
            exp_min_e = params.fixed_columns::<U1>(1)
                              .component_mul(&params.fixed_columns::<U1>(2));
            sincos_phi = Self::random_unit_3x2d(rng);
        } else {
            // Compatibility-oriented version
            cos_theta = Vector3::zero();
            exp_min_e = Vector3::zero();
            sincos_phi = Matrix3x2::zero();
            for par in 0..OUTGOING_COUNT {
                cos_theta[par] = 2. * rng.random() - 1.;
                let phi = 2. * PI * rng.random();
                sincos_phi[(par, X)] = sin(phi);
                sincos_phi[(par, Y)] = cos(phi);
                exp_min_e[par] = rng.random() * rng.random();
            }
        }

        // Compute the momenta
        let sin_theta = cos_theta.map(|cos| sqrt(1. - sqr(cos)));
        let energy = exp_min_e.map(|e_me| -ln(e_me));
        Matrix4x3::from_fn(|coord, par| {
            energy[par] * match coord {
                X => sin_theta[par] * sincos_phi[(par, X)],
                Y => sin_theta[par] * sincos_phi[(par, Y)],
                Z => cos_theta[par],
                E => 1.,
                _ => unreachable!()
            }
        })
    }

    /// Generate vectors on the unit circle with uniform angle distribution
    ///
    /// NOTE: Similar techniques may be used to generate a vector on the unit
    ///       sphere, but that benchmarked unfavorably, likely because it
    ///       entails bringing more computations close to the RNG calls and
    ///       because the 2D case fits available vector hardware more tightly.
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
                    X => sin(phi[par]),
                    Y => cos(phi[par]),
                    _ => unreachable!()
                }
            })
        }
    }

    /// Simulate the impact of N calls to "generate()" on an RNG
    #[cfg(all(feature = "multi-threading",
              not(feature = "faster-threading")))]
    pub(crate) fn simulate_event_batch(rng: &mut RandomGenerator,
                                       num_events: usize) {
        if cfg!(all(feature = "simd-friendly-rng", feature = "fast-sincos")) {
            // If fast-sincos is enabled, the number of RNG calls per event is
            // nondeterministic, so we must simulate events one by one.
            for _ in 0..num_events*OUTGOING_COUNT {
                rng.skip(9);
                Self::random_unit_3x2d(rng);
            }
        } else {
            // Otherwise, we know how many RNG calls will be made per event
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

    /// Extract the 4-momentum of a single particle (internal for now)
    fn momentum(&self, par: usize) -> Momentum {
        Momentum::from_iterator(self.0.fixed_rows::<U1>(par).iter().cloned())
    }

    /// Extract the electron 4-momentum
    pub fn electron_momentum(&self) -> Momentum {
        self.momentum(INCOMING_E_M)
    }

    /// Extract the positron 4-momentum
    #[allow(dead_code)]
    pub fn positron_momentum(&self) -> Momentum {
        self.momentum(INCOMING_E_P)
    }

    /// Extract a photon's 4-momentum
    pub fn outgoing_momentum(&self, par: usize) -> Momentum {
        self.momentum(OUTGOING_SHIFT + par)
    }

    /// Access the outgoing 4-momenta
    pub fn outgoing_momenta(&self) -> OutgoingMomentaSlice {
        self.0.fixed_rows::<U3>(OUTGOING_SHIFT)
    }

    /// Minimal outgoing photon energy
    pub fn min_photon_energy(&self) -> Real {
        if cfg!(feature = "no-photon-sorting") {
            let first_out_e = self.outgoing_momenta()[(0, E)];
            self.outgoing_momenta()
                .fixed_columns::<U1>(E)
                .iter()
                .skip(1)
                .fold(first_out_e, |e1, &e2| if e1 < e2 { e1 } else { e2 })
        } else {
            // Use the fact that photons are sorted by decreasing energy
            self.outgoing_momenta()[(OUTGOING_COUNT-1, E)]
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