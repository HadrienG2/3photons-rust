//! Facilities for computing spinor products

use ::{
    event::{Event, INCOMING_COUNT, OUTGOING_COUNT, PARTICLE_COUNT},
    linalg::{X, Y, Z, E},
    numeric::{
        Complex,
        functions::{sqr, sqrt},
        Real,
        reals::consts::SQRT_2,
    },
};

use num_traits::Zero;


/// Square root of eight
const RAC8: Real = 2. * SQRT_2;


/// Massless 4-momenta spinor inner products
pub struct SpinorProducts {
    /// Gram matrix associated with the inner products
    sx: [[Complex; PARTICLE_COUNT]; PARTICLE_COUNT],
}
//
impl SpinorProducts {
    // ### CONSTRUCTION ###

    /// Build spinor products from previously generated particle 4-momenta
    pub fn new(event: &Event) -> Self {
        // The underlying GramMatrix is not particularly specific to the
        // physical problem of e+e- -> ppp collisions, but this struct is
        debug_assert_eq!(INCOMING_COUNT, 2);
        debug_assert_eq!(OUTGOING_COUNT, 3);

        // Access the array of incoming and outgoing particle 4-momenta
        let p_arr = event.all_momenta();

        // Compute the spinor products (method from M. Mangano and S. Parke)
        // TODO: Once Rust supports it, initialize xx and fx more directly
        let mut xx = [0.; PARTICLE_COUNT];
        let mut fx = [Complex::new(0., 0.); PARTICLE_COUNT];
        for (i, p) in p_arr.iter().enumerate() {
            xx[i] = sqrt(p[E] + p[Z]);
            fx[i] = Complex::new(p[X], p[Y]) / xx[i];
        }

        // Fill up the matrix
        let mut result = Self {
            sx: [[Complex::zero(); PARTICLE_COUNT]; PARTICLE_COUNT],
        };
        for j in 1..PARTICLE_COUNT {
            for i in 0..j {
                let cx = fx[i]*xx[j] - fx[j]*xx[i];
                result.sx[i][j] = cx;
                result.sx[j][i] = -cx;
            }
        }

        // Return the result
        result
    }


    // ### GRAM MATRIX ACCESSORS ###

    /// Quickly access the underlying Gram matrix of spinor products
    #[inline]
    fn s(&self, i: usize, j: usize) -> Complex {
        self.sx[i][j]
    }

    #[inline]
    fn t(&self, i: usize, j: usize) -> Complex {
        -self.s(i, j).conj()
    }

    #[inline]
    pub fn s2(&self, i: usize, j: usize) -> Real {
        self.s(i, j).norm_sqr()
    }


    // ### AMPLITUDE COMPUTATIONS ###

    /// Standard amplitude for helicities ++-
    #[inline]
    pub fn a_ppm(&self, k1: usize, k2: usize, k3: usize) -> Complex {
        -RAC8 * self.s(0, 1) * sqr(self.s(0, k3)) /
            (self.s(0, k1) * self.s(0, k2) * self.s(1, k1) * self.s(1, k2))
    }

    /// Standard amplitude for helicities +--
    #[inline]
    pub fn a_pmm(&self, k1: usize, k2: usize, k3: usize) -> Complex {
        -RAC8 * self.t(0, 1) * sqr(self.t(1, k1)) /
            (self.t(1, k2) * self.t(1, k3) * self.t(0, k2) * self.t(0, k3))
    }

    /// Anomalous amplitude for helicities ++-
    #[inline]
    pub fn b_ppm(&self, k1: usize, k2: usize, k3: usize) -> Complex {
        -RAC8 * self.t(0, 1) * sqr(self.t(k1, k2) * self.s(k3, 0))
    }

    /// Anomalous amplitude for helicities +--
    #[inline]
    pub fn b_pmm(&self, k1: usize, k2: usize, k3: usize) -> Complex {
        -RAC8 * self.s(0, 1) * sqr(self.t(k1, 1) * self.s(k2, k3))
    }

    /// Anomalous amplitude for helicities +++
    #[inline]
    pub fn b_ppp(&self, k1: usize, k2: usize, k3: usize) -> Complex {
        -RAC8 * self.s(0, 1) * (sqr(self.t(k1, k2) * self.t(k3, 1)) +
                                sqr(self.t(k1, k3) * self.t(k2, 1)) +
                                sqr(self.t(k2, k3) * self.t(k1, 1)))
    }

    /// Anomalous amplitude for helicities ---
    #[inline]
    pub fn b_mmm(&self, k1: usize, k2: usize, k3: usize) -> Complex {
        -RAC8 * self.t(0, 1) * (sqr(self.s(k1, 0) * self.s(k2, k3)) +
                                sqr(self.s(k2, 0) * self.s(k1, k3)) +
                                sqr(self.s(k3, 0) * self.s(k1, k2)))
    }
}