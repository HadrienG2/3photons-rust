//! Some shared linear algebra concepts

use ::numeric::Real;

use nalgebra::{
    core::dimension::*,
    MatrixSlice,
    MatrixSliceMut,
    Vector4,
    VectorN,
};


// ### BASIC VECTOR TYPES ###

/// Re-export of some nalgebra types
pub use nalgebra::{Vector2, Vector5, Matrix5};

/// Re-export of nalgebra's 8-vector type
pub type Vector8<T> = VectorN<T, U8>;


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
pub fn xyz(m: &Momentum) -> MatrixSlice<Real, U3, U1, U1, U4> {
    m.fixed_rows::<U3>(X)
}

/// Get a mutable view on the spatial part of a 4-momentum
pub fn xyz_mut(m: &mut Momentum) -> MatrixSliceMut<Real, U3, U1, U1, U4> {
    m.fixed_rows_mut::<U3>(X)
}