//! Facilities for computing spinor products

use ::{
    event::{
        Event,
        INCOMING_COUNT,
        INCOMING_E_M as E_M,
        INCOMING_E_P as E_P,
        OUTGOING_COUNT,
        ParticleMatrix,
    },
    linalg::{E, U1, X, Y, Z},
    numeric::{
        Complex,
        functions::{conj, sqr, sqrt},
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
        let p = event.all_momenta();
        let p_x = p.fixed_columns::<U1>(X);
        let p_y = p.fixed_columns::<U1>(Y);
        let p_z = p.fixed_columns::<U1>(Z);
        let p_e = p.fixed_columns::<U1>(E);

        // Compute the spinor products (method from M. Mangano and S. Parke)
        let xx = (p_e + p_z).map(sqrt);
        let inv_xx = xx.map(|x| 1. / x);
        let re_fx = p_x.component_mul(&inv_xx);
        let im_fx = p_y.component_mul(&inv_xx);

        // Fill up the Gram matrix
        // TODO: Can we leverage antisymmetry + zero diagonal better?
        let result = Self {
            sx: ParticleMatrix::from_fn(|i, j| {
                Complex::new(re_fx[i]*xx[j] - re_fx[j]*xx[i],
                             im_fx[i]*xx[j] - im_fx[j]*xx[i])
            }),
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
        -conj(self.s(i, j))
    }


    // ### AMPLITUDE COMPUTATIONS ###
    //
    // TODO: These long products of complex numbers are getting expensive,
    //       try expanding them on paper to see if they can be reworded into a
    //       simpler overall expression.
    //
    // TODO: Also, the RAC8 prefactor should be applied to the final result, not
    //       to every contribution thereof.

    /// Standard amplitude for given photon helicities
    #[inline]
    pub fn a(&self, helicities: PhotonHelicities) -> Complex {
        use self::PhotonHelicities::*;
        match helicities {
            MMM => Complex::zero(),
            MMP => self.a_pmm(4, 2, 3),
            MPM => self.a_pmm(3, 2, 4),
            MPP => self.a_ppm(3, 4, 2),
            PMM => self.a_pmm(2, 3, 4),
            PMP => self.a_ppm(2, 4, 3),
            PPM => self.a_ppm(2, 3, 4),
            PPP => Complex::zero(),
        }
    }

    /// Standard amplitude for helicities ++-
    #[inline]
    fn a_ppm(&self, k1: usize, k2: usize, k3: usize) -> Complex {
        -RAC8 * self.s(E_M, E_P) * sqr(self.s(E_M, k3)) /
            (self.s(E_M, k1) * self.s(E_M, k2) * self.s(E_P, k1) * self.s(E_P, k2))
    }

    /// Standard amplitude for helicities +--
    #[inline]
    fn a_pmm(&self, k1: usize, k2: usize, k3: usize) -> Complex {
        -RAC8 * self.t(E_M, E_P) * sqr(self.t(E_P, k1)) /
            (self.t(E_P, k2) * self.t(E_P, k3) * self.t(E_M, k2) * self.t(E_M, k3))
    }

    /// Anomalous amplitude 𝛽₊ for given photon helicities
    #[inline]
    pub fn b_p(&self, helicities: PhotonHelicities) -> Complex {
        use self::PhotonHelicities::*;
        match helicities {
            MMM => Complex::zero(),
            MMP => self.bp_pmm(4, 2, 3),
            MPM => self.bp_pmm(3, 2, 4),
            MPP => self.bp_ppm(3, 4, 2),
            PMM => self.bp_pmm(2, 3, 4),
            PMP => self.bp_ppm(2, 4, 3),
            PPM => self.bp_ppm(2, 3, 4),
            PPP => Complex::zero(),
        }
    }

    /// Anomalous amplitude 𝛽₊ for helicities ++-
    #[inline]
    fn bp_ppm(&self, k1: usize, k2: usize, k3: usize) -> Complex {
        -RAC8 * self.t(E_M, E_P) * sqr(self.t(k1, k2) * self.s(k3, E_M))
    }

    /// Anomalous amplitude 𝛽₊ for helicities +--
    #[inline]
    fn bp_pmm(&self, k1: usize, k2: usize, k3: usize) -> Complex {
        -RAC8 * self.s(E_M, E_P) * sqr(self.t(k1, E_P) * self.s(k2, k3))
    }

    /// Anomalous amplitude 𝛽₋ for given photon helicities
    #[inline]
    pub fn b_m(&self, helicities: PhotonHelicities) -> Complex {
        use self::PhotonHelicities::*;
        match helicities {
            MMM => self.bm_mmm(2, 3, 4),
            MMP => Complex::zero(),
            MPM => Complex::zero(),
            MPP => Complex::zero(),
            PMM => Complex::zero(),
            PMP => Complex::zero(),
            PPM => Complex::zero(),
            PPP => self.bm_ppp(2, 3, 4),
        }
    }

    /// Anomalous amplitude 𝛽₋ for helicities +++
    #[inline]
    fn bm_ppp(&self, k1: usize, k2: usize, k3: usize) -> Complex {
        -RAC8 * self.s(E_M, E_P) * (sqr(self.t(k1, k2) * self.t(k3, E_P)) +
                                    sqr(self.t(k1, k3) * self.t(k2, E_P)) +
                                    sqr(self.t(k2, k3) * self.t(k1, E_P)))
    }

    /// Anomalous amplitude 𝛽₋ for helicities ---
    #[inline]
    fn bm_mmm(&self, k1: usize, k2: usize, k3: usize) -> Complex {
        -RAC8 * self.t(E_M, E_P) * (sqr(self.s(k1, E_M) * self.s(k2, k3)) +
                                    sqr(self.s(k2, E_M) * self.s(k1, k3)) +
                                    sqr(self.s(k3, E_M) * self.s(k1, k2)))
    }
}


/// Output photon helicities (M is - and P is +)
#[allow(missing_docs)]
#[derive(Clone, Copy, Debug, Eq, PartialEq)]
pub enum PhotonHelicities {
    MMM,
    MMP,
    MPM,
    MPP,
    PMM,
    PMP,
    PPM,
    PPP,
}