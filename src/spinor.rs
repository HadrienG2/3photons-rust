//! Facilities for computing spinor products

use event::{Event, INCOMING_COUNT, OUTGOING_COUNT, PARTICLE_COUNT};
use linalg::{X, Y, Z, E};
use numeric::{Complex, Real, sqr, sqrt};
use numeric::reals::consts::SQRT_2;
use num_traits::Zero;


/// Square root of eight
const RAC8: Real = 2. * SQRT_2;


/// Array of massless momenta spinor products
pub struct SpinorProducts {
    /// Massless momenta spinor inner products Gram matrix
    gram_matrix: ComplexGramMatrix,
}
//
impl SpinorProducts {
    // ### CONSTRUCTION ###

    /// Build spinor products from previously generated particle momenta
    pub fn new(event: &Event) -> Self {
        // The underlying GramMatrix is not particularly specific to the
        // physical problem of e+e- -> ppp collisions, but this struct is...
        debug_assert_eq!(INCOMING_COUNT, 2);
        debug_assert_eq!(OUTGOING_COUNT, 3);

        // ...even though all it does is computations on a GramMatrix
        SpinorProducts {
            gram_matrix: ComplexGramMatrix::new(event),
        }
    }


    // ### GRAM MATRIX ACCESSORS ###

    /// Quickly access the underlying Gram matrix of spinor products
    fn s(&self, i: usize, j: usize) -> Complex {
        self.gram_matrix.s(i, j)
    }

    fn t(&self, i: usize, j: usize) -> Complex {
        self.gram_matrix.t(i, j)
    }

    pub fn s2(&self, i: usize, j: usize) -> Real {
        self.s(i, j).norm_sqr()
    }


    // ### AMPLITUDE COMPUTATIONS ###

    /// Standard amplitude for helicities ++-
    pub fn a_ppm(&self, k1: usize, k2: usize, k3: usize) -> Complex {
        -RAC8 * self.s(0, 1) * sqr(self.s(0, k3)) /
            (self.s(0, k1) * self.s(0, k2) * self.s(1, k1) * self.s(1, k2))
    }

    /// Standard amplitude for helicities +--
    pub fn a_pmm(&self, k1: usize, k2: usize, k3: usize) -> Complex {
        -RAC8 * self.t(0, 1) * sqr(self.t(1, k1)) /
            (self.t(1, k2) * self.t(1, k3) * self.t(0, k2) * self.t(0, k3))
    }

    /// Anomalous amplitude for helicities ++-
    pub fn b_ppm(&self, k1: usize, k2: usize, k3: usize) -> Complex {
        -RAC8 * self.t(0, 1) * sqr(self.t(k1, k2) * self.s(k3, 0))
    }

    /// Anomalous amplitude for helicities +--
    pub fn b_pmm(&self, k1: usize, k2: usize, k3: usize) -> Complex {
        -RAC8 * self.s(0, 1) * sqr(self.t(k1, 1) * self.s(k2, k3))
    }

    /// Anomalous amplitude for helicities +++
    pub fn b_ppp(&self, k1: usize, k2: usize, k3: usize) -> Complex {
        -RAC8 * self.s(0, 1) * (sqr(self.t(k1, k2) * self.t(k3, 1)) +
                                sqr(self.t(k1, k3) * self.t(k2, 1)) +
                                sqr(self.t(k2, k3) * self.t(k1, 1)))
    }

    /// Anomalous amplitude for helicities ---
    pub fn b_mmm(&self, k1: usize, k2: usize, k3: usize) -> Complex {
        -RAC8 * self.t(0, 1) * (sqr(self.s(k1, 0) * self.s(k2, k3)) +
                                sqr(self.s(k2, 0) * self.s(k1, k3)) +
                                sqr(self.s(k3, 0) * self.s(k1, k2)))
    }
}


// ### IMPLEMENTATION DETAILS ###

/// Gram matrix of massless momenta direct & conjugate spinor inner products.
/// Extracted from SpinorProducts as it uses an optimized but unpleasant data
/// storage model which is best left on its own.
struct ComplexGramMatrix {
    /// You should ONLY access this data through the s(i, j) method. The reason
    /// is that we only fill the lower half of the matrix, and leave the rest of
    /// the matrix full of garbage, opting to compute the remaining values on
    /// the fly if requested instead of storing them.
    ///
    /// We found this to improve performance over the naive approach of filling
    /// and accessing the entire Gram matrix.
    ///
    /// We explored using packed storage instead, but performance analysis found
    /// the extra indexing complexity to be a worse trade-off.
    ///
    /// We also explored computing the whole contents of the matrix on the fly,
    /// and again found this to be a worse trade-off.
    ///
    /// So it looks like this is the most efficient storage model.
    ///
    sx: [[Complex; PARTICLE_COUNT]; PARTICLE_COUNT],
}
//
impl ComplexGramMatrix {
    /// Build a Gram matrix from previously generated particle momenta
    fn new(event: &Event) -> Self {
        // Access the array of incoming and outgoing particle momenta
        let p_arr = event.dump_momenta();

        // Compute the spinor products (method from M. Mangano and S. Parke)
        // TODO: Once Rust supports it, initialize xx and fx more directly
        let mut xx = [0.; PARTICLE_COUNT];
        let mut fx = [Complex::new(0., 0.); PARTICLE_COUNT];
        for (i, p) in p_arr.iter().enumerate() {
            xx[i] = sqrt(p[E] + p[Z]);
            fx[i] = Complex::new(p[X], p[Y]) / xx[i];
        }

        // We will only fill the lower half of the matrix
        let mut result = ComplexGramMatrix {
            sx: [[Complex::zero(); PARTICLE_COUNT]; PARTICLE_COUNT],
        };
        for j in 1..PARTICLE_COUNT {
            for i in 0..j {
                let cx = fx[i]*xx[j] - fx[j]*xx[i];
                result.sx[i][j] = cx;
            }
        }

        // Return the result
        result
    }

    /// Massless momenta spinor inner products Gram matrix accessor
    fn s(&self, i: usize, j: usize) -> Complex {
        if i < j {
            self.sx[i][j]
        } else {
            -self.sx[j][i]
        }
    }

    /// Massless momenta conjugate spinor inner products Gram matrix accessor
    fn t(&self, i: usize, j: usize) -> Complex {
        -self.s(i, j).conj()
    }
}
