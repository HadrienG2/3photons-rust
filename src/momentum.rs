//! This module implements some domain-specific 4-momentum handling logic.

use crate::numeric::Float;
use nalgebra::SVector;

/// 4-momentum dimension
pub const MOMENTUM_DIM: usize = 4;

/// Relativistic 4-momentum
pub type Momentum = SVector<Float, MOMENTUM_DIM>;

/// Convenience const for accessing the X coordinate of a 4-vector
pub const X: usize = 0;

/// Convenience const for accessing the Y coordinate of a 4-vector
pub const Y: usize = 1;

/// Convenience const for accessing the Z coordinate of a 4-vector
pub const Z: usize = 2;

/// Convenience const for accessing the E coordinate of a 4-vector
pub const E: usize = 3;
