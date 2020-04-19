//! Manipulation of matrix elements

use crate::{
    coupling::Couplings,
    event::{Event, NUM_OUTGOING, NUM_SPINS},
    linalg::vecmat::*,
    numeric::{functions::*, Float},
    spinor::SpinorProducts,
};
use std::fmt::Display;

// ### MATRIX ELEMENTS ###

/// Number of matrix elements
pub const NUM_MAT_ELEMS: usize = 5;

/// Storage for per-matrix element data
pub type MEsVector = Vector5<Float>;

/// Index of the electromagnetic element
pub const A: usize = 0;

/// Index of the positive electroweak element
pub const B_P: usize = 1;

/// Index of the negative electroweak element
pub const B_M: usize = 2;

/// Index of the real part of the mixed element
pub const R_MX: usize = 3;

/// Index of the imaginary part of the mixed element
pub const I_MX: usize = 4;

// ### PER-HELICITY CONTRIBUTIONS TO MATRIX ELEMENTS ###

/// Array of square matrix elements contribution with detail of helicities
pub struct MEsContributions {
    /// Array of squared matrix elements, featuring five contributions with the
    /// detail of outgoing helicities configuration
    ///
    /// The rows are the same as in ResultVector, the columns map to spin
    /// configurations encoded as a binary number:
    ///     - Configuration 0 (0b000) is ---
    ///     - Configuration 1 (0b001) is --+
    ///     - And so on...
    ///
    m2: Matrix5x8<Float>,
}
//
impl MEsContributions {
    /// Construct the matrix element contributions from the spinor products
    pub fn new(couplings: &Couplings, event: &Event) -> Self {
        // This code is very specific to the current problem definition
        assert_eq!(NUM_OUTGOING, 3);
        assert_eq!(NUM_MAT_ELEMS, 5);

        // Compute spinor inner products
        let spinor = SpinorProducts::new(&event);

        // Compute the helicity amplitudes, formerly known as a, b_p and b_m,
        // for each possible output spin configuration
        use crate::spinor::PhotonHelicities::*;
        let helicities = Vector8::from_column_slice(&[MMM, MMP, MPM, MPP, PMM, PMP, PPM, PPP]);
        let a_amps = helicities.map(|hel| spinor.a(hel) * couplings.g_a);
        let bp_amps = helicities.map(|hel| spinor.b_p(hel) * couplings.g_beta_p);
        let bm_amps = helicities.map(|hel| spinor.b_m(hel) * couplings.g_beta_m);
        let mixed_amps = a_amps.zip_map(&bp_amps, |a, b_p| 2. * a * conj(b_p));

        // Compute the matrix elements
        MEsContributions {
            m2: Matrix5x8::from_fn(|contrib, hel| match contrib {
                A => norm_sqr(a_amps[hel]),
                B_P => norm_sqr(bp_amps[hel]),
                B_M => norm_sqr(bm_amps[hel]),
                R_MX => re(mixed_amps[hel]),
                I_MX => im(mixed_amps[hel]),
                _ => unreachable!(),
            }),
        }
    }

    /// Compute the sums of the squared matrix elements for each contribution
    pub fn m2_sums(&self) -> MEsVector {
        self.m2.column_sum()
    }
}

impl Display for MEsContributions {
    fn fmt(&self, fmt: &mut std::fmt::Formatter<'_>) -> std::fmt::Result {
        assert_eq!(NUM_OUTGOING, 3);
        assert_eq!(NUM_SPINS, 2);

        for index in 0..NUM_MAT_ELEMS {
            writeln!(fmt, "Matrix element #{}", index)?;
            writeln!(fmt, "---  \t--+  \t-+-  \t-++  \t+--  \t+-+  \t++-  \t+++")?;
            let contribution = self.m2.row(index);
            for &matrix_elem in contribution.iter() {
                write!(fmt, "{}  \t", matrix_elem)?;
            }
            writeln!(fmt)?;
        }
        Ok(())
    }
}
