//! Intermediary results of the computation for one event

use ::{
    coupling::Couplings,
    event::{Event, OUTGOING_COUNT},
    linalg::{Matrix5x8, U1, Vector5, Vector8},
    numeric::{
        functions::conj,
        Real
    },
    spinor::SpinorProducts,
};


/// Number of results (matrix elements)
pub const NUM_RESULTS: usize = 5;

/// Storage for matrix elements
pub type ResultVector<T> = Vector5<T>;

/// Index of the electromagnetic matrix element
pub const A: usize = 0;

/// Index of the positive electroweak matrix element
pub const B_P: usize = 1;

/// Index of the negative electroweak matrix element
pub const B_M: usize = 2;

/// Index of the real part of the mixed matrix element
pub const R_MX: usize = 3;

/// Index of the imaginary part of the mixed matrix element
pub const I_MX: usize = 4;


/// Array of square matrix elements contribution with detail of helicities
pub struct ResultContribution {
    /// Array of squared matrix elements, featuring five contributions with the
    /// detail of outgoing helicities configuration
    ///
    /// The rows are the same as in ResultVector, the columns map to spin
    /// configurations encoded as a binary number:
    ///     - Configuration 0 (0b000) is ---
    ///     - Configuration 1 (0b001) is --+
    ///     - And so on...
    ///
    m2x: Matrix5x8<Real>,
}
//
impl ResultContribution {
    /// Construct the matrix element from the spinor products
    pub fn new(couplings: &Couplings, event: Event) -> Self {
        // This code is very specific to the current problem definition
        assert_eq!(OUTGOING_COUNT, 3);
        assert_eq!(NUM_RESULTS, 5);

        // Compute spinor inner products
        let spinor = SpinorProducts::new(event);

        // Compute the helicity amplitudes, formerly known as a, b_p and b_m,
        // for each possible output spin configuration
        //
        // TODO: These zeroes might be wasted computations, try doing without
        //       them and see if it helps.
        //
        // TODO: Review data layout (e.g. tuples vs array, row- vs col-major...)
        //
        use spinor::PhotonHelicities::*;
        let helicities = Vector8::from_column_slice(
            &[MMM, MMP, MPM, MPP, PMM, PMP, PPM, PPP]
        );
        let a_amps = helicities.map(|hel| spinor.a(hel) * couplings.g_a);
        let bp_amps = helicities.map(|hel| spinor.b_p(hel) * couplings.g_bp);
        let bm_amps = helicities.map(|hel| spinor.b_m(hel) * couplings.g_bm);
        let mixed_amp = Vector8::from_fn(|hel, _| {
            2. * a_amps[hel] * conj(bp_amps[hel])
        });

        // Compute the matrix elements
        ResultContribution {
            m2x: Matrix5x8::from_fn(|contrib, hel| {
                match contrib {
                    A => a_amps[hel].norm_sqr(),
                    B_P => bp_amps[hel].norm_sqr(),
                    B_M => bm_amps[hel].norm_sqr(),
                    R_MX => mixed_amp[hel].re,
                    I_MX => mixed_amp[hel].im,
                    _ => unreachable!(),
                }
            }),
        }
    }

    /// Compute the sums of the squared matrix elements for each contribution
    pub fn m2_sums(&self) -> ResultVector<Real> {
        ResultVector::from_fn(|i, _| self.m2x.fixed_rows::<U1>(i).iter().sum())
    }

    /// Display the results in human-readable form
    #[allow(dead_code)]
    pub fn display(&self) {
        assert_eq!(OUTGOING_COUNT, 3);
        
        for index in 0..NUM_RESULTS {
            println!("Contribution {}", index);
            println!("---  \t--+  \t-+-  \t-++  \t+--  \t+-+  \t++-  \t+++");
            let contribution = self.m2x.fixed_rows::<U1>(index);
            for &matrix_elem in contribution.iter() {
                print!("{}  \t", matrix_elem);
            }
            println!();
        }
    }
}
