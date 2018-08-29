//! This module takes care of event generation and storage

use ::{
    linalg::{self, Momentum, X, Y, Z, E},
    numeric::{
        cos,
        exp,
        ln,
        Real,
        reals::consts::{PI, FRAC_PI_2},
        sin,
        sqr,
        sqrt
    },
    ranf::RanfGenerator,
};

use num_traits::identities::Zero;


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
pub struct Event {
    /// Array of incoming and outgoing 4-momenta
    p: Momenta,

    /// Total center-of-mass energy of the collision
    e_tot: Real,

    /// Weight of generated events
    ev_weight: Real,
}
//
impl Event {
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
        let weight = (2. * (OUTGOING_COUNT as Real) - 4.) * ln(e_tot) + z;
        assert!((weight >= -180.) && (weight <= 174.));
        let weight = exp(weight);
        
        // Construct and return the output data structure
        let mut result = Event {
            p: [Momentum::zero(); PARTICLE_COUNT],
            e_tot: e_tot,
            ev_weight: weight,
        };
        result.p[INCOMING_E_M] = Momentum::new(-e_tot/2., 0., 0., e_tot/2.);
        result.p[INCOMING_E_P] = Momentum::new(e_tot/2., 0., 0., e_tot/2.);
        // NOTE: The outgoing momenta have been zero-initialized above
        result
    }


    // ### EVENT GENERATION ###

    /// Use a highly specialized version of the RAMBO (RAndom Momenta
    /// Beautifully Organized) algorithm from S.D. Ellis, R. Kleiss and W.J.
    /// Stirling in order to generate the momenta of the three outgoing photons.
    ///
    /// All events have the same weight, which can be read via self.ev_weight().
    ///
    /// "rng" is a random number generator. For now, it has to be RANF, later
    /// we'll make the code generic with respect to this consideration.
    ///
    pub fn generate_momenta(&mut self, rng: &mut RanfGenerator) {
        // Generate massless outgoing momenta in infinite phase space
        // TODO: Once Rust supports it, initialize q_arr more directly
        let mut q_arr = [Momentum::zero(); OUTGOING_COUNT];
        for q in q_arr.iter_mut() {
            let cos_theta = 2. * rng.random() - 1.;
            let sin_theta = sqrt(1.0 - sqr(cos_theta));
            let phi = 2. * PI * rng.random();
            q[E] = - ln(rng.random() * rng.random());
            q[Z] = q[E] * cos_theta;
            q[Y] = q[E] * sin_theta * cos(phi);
            q[X] = q[E] * sin_theta * sin(phi);
        }
        let q_arr = q_arr;

        // Calculate the parameters of the conformal transformation
        let r = q_arr.iter().sum();
        let r_norm = 1. / linalg::lorentz_norm(&r);
        let b = -linalg::xyz(&r) * r_norm;
        let g = r[E] * r_norm;
        let a = 1. / (1. + g);
        let d = self.e_tot * r_norm;  // NOTE: Was called X in original code

        // Transform the Q's conformally into output momenta
        let p_out = self.dump_outgoing_mut();
        for (p, q) in p_out.iter_mut().zip(q_arr.iter()) {
            let q_xyz = linalg::xyz(q);
            let bq = b.dot(&q_xyz);
            let p_xyz = d * (q_xyz + b * (q[E] + a * bq));
            let p_e = d * (g * q[E] + bq);
            *p = Momentum::new(p_xyz[X], p_xyz[Y], p_xyz[Z], p_e);
        }
    }


    // ### PHOTON SORTING ###

    /// Sort outgoing photons in order of decreasing energy
    /// Roughly equivalent to the original ppp::TRI method
    pub fn sort_output_momenta(&mut self) {
        // As currently written, this function only works with 3 output momenta
        debug_assert_eq!(OUTGOING_COUNT, 3);

        // Original sorting algorithm from the C++ version of 3photons
        self.sort_output_indices(0, 1);
        self.sort_output_indices(0, 2);
        self.sort_output_indices(1, 2);
    }

    /// Order output momenta i and j, where it is assumed that i < j, according
    /// to their energy. This is just a helper for the sorting method.
    fn sort_output_indices(&mut self, i: usize, j: usize) {
        debug_assert!(i < j);
        let p_out = self.dump_outgoing_mut();
        if p_out[j][E] > p_out[i][E] {
            p_out.swap(i, j);
        }
    }


    // ### ACCESSORS ###

    /// Access the event weight
    pub fn weight(&self) -> Real {
        self.ev_weight
    }

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
