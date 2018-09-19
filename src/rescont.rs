//! Intermediary results of the computation for one event

use ::{
    coupling::Couplings,
    event::{Event, INCOMING_COUNT, OUTGOING_COUNT},
    linalg::{Vector5, Vector8},
    numeric::{Complex, Real},
    spinor::SpinorProducts,
};

use std::mem;


/// Contributions to the squared matrix elements
///
/// The indices of this array map to binary-encoded spin configurations:
///     - Configuration 0 (0b000) is ---
///     - Configuration 1 (0b001) is --+
///     - And so on...
///
type MEContributions = Vector8<Real>;

/// Number of results, and associated storage type
pub const NUM_RESULTS: usize = 5;
pub type ResultVector<T> = Vector5<T>;


/// Array of square matrix elements contribution with detail of helicities
pub struct ResultContribution {
    /// Array of squared matrix elements, featuring five contributions with the
    /// detail of outgoing helicities configuration
    ///
    /// TODO: Should probably be a matrix
    ///
    m2x: ResultVector<MEContributions>,
}
//
impl ResultContribution {
    /// Construct the matrix element from the spinor products
    pub fn new(couplings: &Couplings, event: Event) -> Self {
        // This code is very specific to the current problem definition
        debug_assert_eq!(INCOMING_COUNT, 2);
        debug_assert_eq!(OUTGOING_COUNT, 3);
        debug_assert_eq!(NUM_RESULTS, 5);

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
        let zero = Complex::new(0., 0.);
        let helicity_amps = [
            // ---
            (zero,                 zero,                 spinor.b_mmm(2,3,4)),
            // --+
            (spinor.a_pmm(4,2,3),  spinor.b_pmm(4,2,3),  zero),
            // -+-
            (spinor.a_pmm(3,2,4),  spinor.b_pmm(3,2,4),  zero),
            // -++
            (spinor.a_ppm(3,4,2),  spinor.b_ppm(3,4,2),  zero),
            // +--
            (spinor.a_pmm(2,3,4),  spinor.b_pmm(2,3,4),  zero),
            // +-+
            (spinor.a_ppm(4,2,3),  spinor.b_ppm(4,2,3),  zero),
            // ++-
            (spinor.a_ppm(2,3,4),  spinor.b_ppm(2,3,4),  zero),
            // +++
            (zero,                 zero,                 spinor.b_ppp(2,3,4)),
        ];

        // Compute the matrix elements
        let mut result = ResultContribution {
            m2x: unsafe{ mem::uninitialized() },
        };
        for (index, &(ha, hb_p, hb_m)) in helicity_amps.iter().enumerate() {
            // Take couplings into account
            //
            // TODO: Instead of multiplying by couplings every contribution, try
            //       only applying the couplings to the final result
            //
            let a = ha * couplings.g_a;
            let b_p = hb_p * couplings.g_bp;
            let b_m = hb_m * couplings.g_bm;

            // Compute the squared matrix element terms
            //
            // TODO: This matrix write pattern is inefficient, but maybe that's
            //       the price to pay for an efficient final m2_sums reduction?
            //
            let mixed = 2. * a * b_p.conj();
            result.m2x[0][index] = a.norm_sqr();
            result.m2x[1][index] = b_p.norm_sqr();
            result.m2x[2][index] = b_m.norm_sqr();
            result.m2x[3][index] = mixed.re;
            result.m2x[4][index] = mixed.im;
        }

        // Output the results
        result
    }

    /// Compute the sums of the squared matrix elements for each contribution
    pub fn m2_sums(&self) -> ResultVector<Real> {
        self.m2x.map(|contrib| contrib.iter().sum())
    }

    /// Display the results in human-readable form
    #[allow(dead_code)]
    pub fn display(&self) {
        assert_eq!(OUTGOING_COUNT, 3);
        
        for (index, contribution) in self.m2x.iter().enumerate() {
            println!("Contribution {}", index);
            println!("---  \t--+  \t-+-  \t-++  \t+--  \t+-+  \t++-  \t+++");
            for &matrix_elem in contribution.iter() {
                print!("{}  \t", matrix_elem);
            }
            println!();
        }
    }
}
