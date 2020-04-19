//! Basic numerical concepts used throughout the program

#![allow(missing_docs)]

use num_complex;

// Floating-point precision is configured here
#[cfg(feature = "f32")]
pub type Float = f32;
#[cfg(feature = "f32")]
pub use std::f32 as floats;
#[cfg(not(feature = "f32"))]
pub type Float = f64;
#[cfg(not(feature = "f32"))]
pub use std::f64 as floats;
pub type Complex = num_complex::Complex<Float>;

/// Mathematical functions
pub mod functions {
    /// Compute the conjugate of a Complex number
    pub fn conj(z: super::Complex) -> super::Complex {
        z.conj()
    }

    /// Compute the squared norm of a Complex number
    pub fn norm_sqr(z: super::Complex) -> super::Float {
        z.norm_sqr()
    }

    /// Get the real part of of a Complex number
    pub fn re(z: super::Complex) -> super::Float {
        z.re
    }

    /// Get the imaginary part of of a Complex number
    pub fn im(z: super::Complex) -> super::Float {
        z.im
    }
}
