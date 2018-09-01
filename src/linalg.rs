//! Some shared linear algebra concepts

use ::numeric::Real;

use nalgebra::{
    core::dimension::*,
    self,
    storage::Storage
};


// ### BASIC VECTOR TYPES ###

/// Re-export of nalgebra's 5-vector type
pub use nalgebra::Vector5;

/// Re-export of nalgebra's 8-vector type
pub type Vector8<T> = nalgebra::VectorN<T, U8>;


// ### RELATIVISTIC 4-MOMENTA ###

/// 4-vectors of real numbers, as used by special relativity
type Vector4R = nalgebra::Vector4<Real>;

/// Underlying storage of nalgebra's 4-vector implementation
type V4RImpl = nalgebra::MatrixArray<Real, U4, U1>;

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
pub fn xyz(m: &Momentum)
  -> nalgebra::MatrixSlice<Real, U3, U1,
                           <V4RImpl as Storage<Real, U4, U1>>::RStride,
                           <V4RImpl as Storage<Real, U4, U1>>::CStride>
{
  m.fixed_rows::<U3>(X)
}

/// Get a mutable view on the spatial part of a 4-momentum
pub fn xyz_mut(m: &mut Momentum)
  -> nalgebra::MatrixSliceMut<Real, U3, U1,
                              <V4RImpl as Storage<Real, U4, U1>>::RStride,
                              <V4RImpl as Storage<Real, U4, U1>>::CStride>
{
  m.fixed_rows_mut::<U3>(X)
}