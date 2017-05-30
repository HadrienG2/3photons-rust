//! Facilities for passing around Lorentz scalar products of momenta

use numeric::Real;
use spinor::SpinorProducts;


/// Array of all possible Lorentz 4-momenta products (Gram matrix)
pub struct ScalarProducts<'a> {
    // This matrix is actually so rarely accessed that it's better to compute
    // its elements on the fly from the underlying spinor product matrix
    spinor: &'a SpinorProducts,
}
//
impl<'a> ScalarProducts<'a> {
    /// Build Lorentz scalar products from spinor products
    pub fn new(spinor: &'a SpinorProducts) -> Self {
        ScalarProducts {
            spinor: spinor,
        }
    }

    /// Compute scalar products on the fly
    pub fn ps(&self, i: usize, j: usize) -> Real {
        self.spinor.s2(i, j) / 2.
    }
}
