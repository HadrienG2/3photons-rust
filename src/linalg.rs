//! Some shared linear algebra concepts

use ::numeric::Real;

use nalgebra::{
    core::dimension::*,
    self,
    storage::Storage
};


// ### BASIC VECTOR TYPES ###

/// Re-export some base vector types for further use
pub use nalgebra::Vector5;
pub type Vector8<T> = nalgebra::VectorN<T, U8>;


// ### RELATIVISTIC MOMENTA ###

/// We'll be operating in the space of relativistic 4-momenta
type Vector4 = nalgebra::Vector4<Real>;
type V4Impl = nalgebra::MatrixArray<Real, U4, U1>;
pub type Momentum = Vector4;

/// When manipulating 4-momenta, it may be useful to have more explicit names
/// for the various coordinates.
pub const X: usize = 0;
pub const Y: usize = 1;
pub const Z: usize = 2;
pub const E: usize = 3;

/// Extract the spatial part of a 4-momentum
pub fn xyz(m: &Momentum)
  -> nalgebra::MatrixSlice<Real, U3, U1,
                           <V4Impl as Storage<Real, U4, U1>>::RStride,
                           <V4Impl as Storage<Real, U4, U1>>::CStride>
{
  m.fixed_rows::<U3>(X)
}
pub fn xyz_mut(m: &mut Momentum)
  -> nalgebra::MatrixSliceMut<Real, U3, U1,
                              <V4Impl as Storage<Real, U4, U1>>::RStride,
                              <V4Impl as Storage<Real, U4, U1>>::CStride>
{
  m.fixed_rows_mut::<U3>(X)
}