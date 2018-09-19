//! Facilities for computing spinor products

use ::{
    event::{
        Event,
        INCOMING_COUNT,
        OUTGOING_COUNT,
        ParticleVector,
        ParticleMatrix,
    },
    linalg::{X, Y, Z, E},
    numeric::{
        Complex,
        functions::{sqr, sqrt},
        Real,
        reals::consts::SQRT_2,
    },
};


/// Square root of eight
const RAC8: Real = 2. * SQRT_2;


/// Massless 4-momenta spinor inner products
pub struct SpinorProducts {
    /// Gram matrix associated with the inner products
    sx: ParticleMatrix<Complex>,
}
//
impl SpinorProducts {
    // ### CONSTRUCTION ###

    /// Build spinor products from previously generated particle 4-momenta
    pub fn new(event: Event) -> Self {
        // The underlying Gram matrix is not specific to the physics of
        // e+e- -> ppp collisions, but our methods are specific to it.
        debug_assert_eq!(INCOMING_COUNT, 2);
        debug_assert_eq!(OUTGOING_COUNT, 3);

        // Access the array of incoming and outgoing particle 4-momenta
        let px = event.all_momenta();

        // Compute the spinor products (method from M. Mangano and S. Parke)
        let xx = px.map(|p| sqrt(p[E] + p[Z]));
        let fx = ParticleVector::from_iterator(
            px.iter().zip(xx.iter())
                     .map(|(p, x)| Complex::new(p[X], p[Y]) / x)
        );

        // Fill up the Gram matrix
        // TODO: Can we leverage antisymmetry + zero diagonal better?
        let result = Self {
            sx: ParticleMatrix::from_fn(|i, j| { fx[i]*xx[j] - fx[j]*xx[i] }),
        };

        // Return the result
        result
    }


    // ### GRAM MATRIX ACCESSORS ###

    /// Quickly access the underlying Gram matrix of spinor products
    #[inline]
    fn s(&self, i: usize, j: usize) -> Complex {
        self.sx[(i, j)]
    }

    #[inline]
    fn t(&self, i: usize, j: usize) -> Complex {
        -self.s(i, j).conj()
    }


    // ### AMPLITUDE COMPUTATIONS ###
    //
    // TODO: These long products of complex numbers are getting expensive,
    //       try expanding them on paper to see if they can be reworded into a
    //       simpler overall expression.
    //
    // TODO: Also, the RAC8 prefactor should be applied to the final result, not
    //       to every contribution thereof.

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