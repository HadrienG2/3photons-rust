//! Some shared linear algebra concepts

use ::numeric::Real;

use nalgebra::{
    MatrixSlice,
    MatrixSliceMut,
    Vector4,
    VectorN,
};


// ### BASIC VECTOR TYPES ###

/// Re-export of nalgebra's type-level integers
pub use nalgebra::core::dimension::*;

/// Re-export of some nalgebra types
pub use nalgebra::{Vector2, Vector3, Vector5, Matrix5};

/// Re-export of nalgebra's 8-vector type
pub type Vector8<T> = VectorN<T, U8>;

/// Convenience shorthand for defining vector slices
pub type VectorSlice<'a, T, SliceDim, ParentDim> =
    MatrixSlice<'a, T, SliceDim, U1, U1, ParentDim>;

/// Convenience shorthand for defining mutable vector slices
pub type VectorSliceMut<'a, T, SliceDim, ParentDim> =
    MatrixSliceMut<'a, T, SliceDim, U1, U1, ParentDim>;


// ### RELATIVISTIC 4-MOMENTA ###

/// 4-vectors of real numbers, as used by special relativity
type Vector4R = Vector4<Real>;

/// Relativistic 4-momentum
pub type Momentum = Vector4R;

/// Convenience const for accessing the X coordinate of a 4-vector
pub const X: usize = 0;

/// Convenience const for accessing the Y coordinate of a 4-vector
pub const Y: usize = 1;

/// Convenience const for accessing the Z coordinate of a 4-vector
pub const Z: usize = 2;

/// Convenience const for accessing the E coordinate of a 4-vector
pub const E: usize = 3;

/// Get a read-only view on the spatial part of a 4-momentum
pub fn xyz(m: &Momentum) -> VectorSlice<Real, U3, U4> {
    m.fixed_rows::<U3>(X)
}