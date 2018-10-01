//! Intermediary results of the computation for one event

use ::{
    coupling::Couplings,
    event::{Event, OUTGOING_COUNT},
    linalg::{Matrix5x8, U1, Vector5},
    numeric::{
        functions::conj,
        Real
    },
    spinor::SpinorProducts,
};

use std::mem;


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
        let helicity_amps = [
            [spinor.a(MMM), spinor.b_p(MMM), spinor.b_m(MMM)],
            [spinor.a(MMP), spinor.b_p(MMP), spinor.b_m(MMP)],
            [spinor.a(MPM), spinor.b_p(MPM), spinor.b_m(MPM)],
            [spinor.a(MPP), spinor.b_p(MPP), spinor.b_m(MPP)],
            [spinor.a(PMM), spinor.b_p(PMM), spinor.b_m(PMM)],
            [spinor.a(PMP), spinor.b_p(PMP), spinor.b_m(PMP)],
            [spinor.a(PPM), spinor.b_p(PPM), spinor.b_m(PPM)],
            [spinor.a(PPP), spinor.b_p(PPP), spinor.b_m(PPP)],
        ];

        // Compute the matrix elements
        let mut result = ResultContribution {
            m2x: unsafe{ mem::uninitialized() },
        };
        for (index, &[ha, hb_p, hb_m]) in helicity_amps.iter().enumerate() {
            // Take couplings into account
            let a = ha * couplings.g_a;
            let b_p = hb_p * couplings.g_bp;
            let b_m = hb_m * couplings.g_bm;

            // Compute the squared matrix element terms
            let mixed = 2. * a * conj(b_p);
            result.m2x[(A, index)] = a.norm_sqr();
            result.m2x[(B_P, index)] = b_p.norm_sqr();
            result.m2x[(B_M, index)] = b_m.norm_sqr();
            result.m2x[(R_MX, index)] = mixed.re;
            result.m2x[(I_MX, index)] = mixed.im;
        }

        // Output the results
        result
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
