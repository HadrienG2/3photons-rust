//! Basic numerical concepts used throughout the program

#![allow(missing_docs)]

use num_complex;

// Floating-point precision is configured here
#[cfg(feature = "f32")]
pub type Float = f32;
#[cfg(feature = "f32")]
pub use std::f32 as reals;
#[cfg(not(feature = "f32"))]
pub type Float = f64;
#[cfg(not(feature = "f32"))]
pub use std::f64 as reals;
pub type Complex = num_complex::Complex<Float>;

/// Mathematical functions
pub mod functions {
    /// Compute the conjugate of a Complex number
    pub fn conj(z: super::Complex) -> super::Complex {
        z.conj()
    }
}
