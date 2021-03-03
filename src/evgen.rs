//! This module provides event generation facilities

use crate::{
    event::{Event, NUM_INCOMING, NUM_OUTGOING, NUM_PARTICLES},
    linalg::{
        dimension::*,
        momentum::{E, X, Y, Z},
        vecmat::*,
    },
    numeric::{
        floats::{
            consts::{FRAC_PI_2, PI},
            MIN_POSITIVE,
        },
        Float,
    },
    random::RandomGenerator,
};

use prefix_num_ops::real::*;

/// Generator of ee -> ppp events
pub struct EventGenerator {
    /// Total center-of-mass energy of the collision
    e_total: Float,

    /// Weight of generated events
    ev_weight: Float,

    /// Incoming electron and positron momenta
    incoming_momenta: Matrix2x4<Float>,
}
//
impl EventGenerator {
    // ### CONSTRUCTION ###

    /// Initialize event generation for a center-of-mass energy of e_total.
    ///
    /// Combines former functionality of ppp constructor and IBEGIN-based lazy
    /// initialization from the original C++ 3photons code.
    ///
    #[rustfmt::skip]
    #[allow(clippy::assertions_on_constants)]
    pub fn new(e_total: Float) -> Self {
        // Check on the number of particles. The check for N<101 is gone since
        // unlike the original RAMBO, we don't use arrays of hardcoded size.
        assert!(NUM_OUTGOING > 1);

        // As currently written, this code only works for two incoming particles
        assert_eq!(NUM_INCOMING, 2);

        // Factorials for the phase space weight. Replaces the lazy
        // initialization from the original RAMBO code with less branchy code.
        println!("IBegin");
        // Replaces Z[INP-1] in the original 3photons code
        let mut z_n = ((NUM_OUTGOING - 1) as Float) * ln(FRAC_PI_2);
        for k in 2..NUM_OUTGOING {
            z_n -= 2. * ln((k - 1) as Float);
        }
        let z_n = z_n - ln((NUM_OUTGOING - 1) as Float);

        // NOTE: The check on total energy is gone, because we only generate
        //       massless photons and so the total energy will always be enough.
        //       Counting of nonzero masses is also gone because it was unused.

        // All generated events will have the same weight: pre-compute it
        let ln_weight = (2. * (NUM_OUTGOING as Float) - 4.) * ln(e_total) + z_n;
        assert!((-180. ..=174.).contains(&ln_weight));
        let ev_weight = exp(ln_weight);

        // Compute the incoming particle momenta
        let incoming_momenta = Matrix2x4::new(
            -e_total / 2., 0., 0., e_total / 2.,
            e_total / 2.,  0., 0., e_total / 2.,
        );

        // Construct and return the output data structure
        EventGenerator {
            e_total,
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
        let r = q.column_sum();
        let r_norm_2 = r[E] * r[E] - r.xyz().norm_squared();
        let alpha = self.e_total / r_norm_2;
        let r_norm = sqrt(r_norm_2);
        let beta = 1. / (r_norm + r[E]);

        // Perform the conformal transformation from Q's to output 4-momenta
        let tr_q = q.transpose();
        let tr_q_xyz = tr_q.fixed_columns::<U3>(X);
        let rq = tr_q_xyz * r.xyz();
        let mut p_e = alpha * (r[E] * tr_q.column(E) - rq);
        let b_rq_e = beta * rq - tr_q.column(E);
        let mut p_xyz = alpha * (r_norm * tr_q_xyz + b_rq_e * r.xyz().transpose());

        // Sort the output 4-momenta in order of decreasing energy (if enabled)
        if cfg!(not(feature = "no-photon-sorting")) {
            for par1 in 0..NUM_OUTGOING - 1 {
                for par2 in par1 + 1..NUM_OUTGOING {
                    if p_e[par2] > p_e[par1] {
                        p_e.swap_rows(par1, par2);
                        p_xyz.swap_rows(par1, par2);
                    }
                }
            }
        }

        // Build the final event: incoming momenta + output 4-momenta
        assert_eq!(NUM_PARTICLES, 5, "This part assumes 5-particles events");
        Event::new(Matrix5x4::from_fn(|par, coord| {
            if par < NUM_INCOMING {
                self.incoming_momenta[(par, coord)]
            } else if coord <= Z {
                p_xyz[(par - NUM_INCOMING, coord)]
            } else if coord == E {
                p_e[par - NUM_INCOMING]
            } else {
                unreachable!()
            }
        }))
    }

    /// Generate massless outgoing 4-momenta in infinite phase space
    ///
    /// The output momenta are provided as a matrix where rows are 4-momentum
    /// components (Px, Py, Pz, E) and columns are particles.
    ///
    fn generate_raw(rng: &mut RandomGenerator) -> Matrix4x3<Float> {
        assert_eq!(NUM_OUTGOING, 3, "This part assumes 3 outgoing particles");

        // In all operating modes, random number generation is kept
        // well-separated from computations, as it was observed that it has a
        // harmful interaction with the compiler's loop optimizations.
        if cfg!(feature = "faster-evgen") {
            // This mode allows random number generation to be carried out in a
            // different order, and using different algorithms than what the
            // original 3photons did. This enables greater performance.

            // Generate the basic random parameters of the particles
            let params = Matrix3::from_column_slice(&rng.random9()[..]);
            let cos_theta = params.column(0).map(|r| 2. * r - 1.);
            let exp_min_e = params.column(1).component_mul(&params.column(2));
            let sincos_phi = Self::random_unit_2d_outgoing(rng);

            // Compute the outgoing momenta
            //
            // FIXME: The main obvious remaining bottleneck of this version is
            //        that it spends ~40% of its time computing scalar
            //        logarithms. Using a vectorized ln() implementation in the
            //        computation of the energy vector should help.
            //
            let sin_theta = cos_theta.map(|cos| sqrt(1. - cos.powi(2)));
            let energy = exp_min_e.map(|e_me| -ln(e_me + MIN_POSITIVE));
            Matrix4x3::from_fn(|coord, par| {
                energy[par]
                    * match coord {
                        X => sin_theta[par] * sincos_phi[(par, X)],
                        Y => sin_theta[par] * sincos_phi[(par, Y)],
                        Z => cos_theta[par],
                        E => 1.,
                        _ => unreachable!(),
                    }
            })
        } else {
            // This mode targets maximal reproducibility with respect to the
            // original 3photons program, at the expense of performance.

            // Generate the basic random parameters of the particles
            const COS_THETA: usize = 0;
            const PHI: usize = 1;
            const EXP_MIN_E: usize = 2;
            let params = Matrix3::from_fn(|coord, _par| match coord {
                COS_THETA => 2. * rng.random() - 1.,
                PHI => 2. * PI * rng.random(),
                EXP_MIN_E => rng.random() * rng.random(),
                _ => unreachable!(),
            });
            let cos_theta = params.row(COS_THETA);
            let phi = params.row(PHI);
            let exp_min_e = params.row(EXP_MIN_E);

            // Compute the outgoing momenta
            let cos_phi = phi.map(cos);
            let sin_phi = phi.map(sin);
            let sin_theta = cos_theta.map(|cos| sqrt(1. - cos.powi(2)));
            let energy = exp_min_e.map(|e_me| -ln(e_me + MIN_POSITIVE));
            Matrix4x3::from_fn(|coord, par| {
                energy[par]
                    * match coord {
                        X => sin_theta[par] * sin_phi[par],
                        Y => sin_theta[par] * cos_phi[par],
                        Z => cos_theta[par],
                        E => 1.,
                        _ => unreachable!(),
                    }
            })
        }
    }

    /// Generate 3 vectors on the unit circle with uniform angle distribution
    ///
    /// NOTE: Similar techniques may be used to generate a vector on the unit
    ///       sphere, but that benchmarked unfavorably, likely because...
    ///
    ///       - We perform best in SSE, in which 2 doubles fit better than 3
    ///       - The phi trig ops are expensive, the random cos of theta isn't
    ///       - RNG calls disturb compiler optimizations, and the 3D case brings
    ///         more computations close to them.
    ///       - Statistics force us to discard more points and call the RNG more
    ///
    fn random_unit_2d_outgoing(rng: &mut RandomGenerator) -> Matrix3x2<Float> {
        assert_eq!(NUM_OUTGOING, 3, "This part assumes 3 outgoing particles");

        // Grab three random points on the unit square
        let mut points = Matrix3x2::from_iterator(rng.random6().iter().map(|r| 2. * r - 1.));

        // Re-roll each point until it falls on the unit disc, and is not
        // too close to the origin (otherwise we'll get floating-point issues)
        let mut radii2 = Vector3::from_iterator(points.row_iter().map(|row| row.norm_squared()));
        for (point_idx, radius2) in radii2.iter_mut().enumerate() {
            const MIN_POSITIVE_2: Float = MIN_POSITIVE * MIN_POSITIVE;
            while *radius2 > 1. || *radius2 < MIN_POSITIVE_2 {
                let new_point = Vector2::from_iterator(rng.random2().iter().map(|r| 2. * r - 1.));
                points.set_row(point_idx, &new_point.transpose());
                *radius2 = new_point.norm_squared();
            }
        }

        // Now you only need to normalize to get points on the unit circle
        let norms = radii2.map(|r2| 1. / sqrt(r2));
        for (mut point, norm) in points.row_iter_mut().zip(norms.iter()) {
            point.apply(|coord| coord * norm);
        }
        points
    }

    /// Simulate the impact of N calls to "generate()" on an RNG
    ///
    /// This function must be kept in sync with the `genrate_raw()`
    /// implementation. Such is the price to pay for perfect reproducibility
    /// between single-threaded and multi-threaded runs...
    ///
    #[cfg(all(feature = "multi-threading", not(feature = "faster-threading")))]
    pub fn simulate_event_batch(rng: &mut RandomGenerator, num_events: usize) {
        if cfg!(feature = "faster-evgen") {
            for _ in 0..num_events {
                rng.skip9();
                Self::random_unit_2d_outgoing(rng);
            }
        } else {
            rng.skip(num_events * NUM_OUTGOING * 4);
        }
    }

    // ### EVENT PROPERTIES ###

    /// Access the event weight (identical for all generated events)
    pub fn event_weight(&self) -> Float {
        self.ev_weight
    }
}
