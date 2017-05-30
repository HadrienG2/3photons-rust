//! Basic numerical concepts used throughout the program

use num_complex;
use num_traits::One;
use std::ops::{Div, Mul};


/// Floating-point precision is configured here
pub use std::f64 as reals;
pub type Real = f64;
pub type Complex = num_complex::Complex<Real>;


// The following functions are trivial wrappers around member functions of
// floating-point types. Their only purpose is to achieve a consistent notation
// in floating-point expressions, instead of the unpleasant mixture of standard
// calculator notation and RPN that the Rust standard library forces us into.

/// Absolute value of a real number
pub fn abs(x: Real) -> Real {
    x.abs()
}

/// Square root of a real number
pub fn sqrt(x: Real) -> Real {
    x.sqrt()
}

/// Natural logarithm of a real number
pub fn ln(x: Real) -> Real {
    x.ln()
}

/// Exponential of a real number
pub fn exp(x: Real) -> Real {
    x.exp()
}

/// Sine of a real number
pub fn sin(x: Real) -> Real {
    x.sin()
}

/// Cosine of a real number
pub fn cos(x: Real) -> Real {
    x.cos()
}


// The following functions are generic in nature, and designed to work on both
// real and complex numbers.

/// Compute the square of any number, optimized shorthand for powi(x, 2)
pub fn sqr<T>(x: T) -> T
    where T: Mul<Output=T> + Copy
{
    x * x
}

/// Raise a number to an arbitrary integer power, like {float}.powi(), using a
/// binary exponentiation algorithm.
///
/// This recursive version is clean and concise, but should only be used when n
/// is known at compile time, so that the compiler can unroll the recursion...
///
pub fn powi<T>(x: T, n: i32) -> T
    where T: Mul<Output=T> + Div<Output=T> + Copy + One
{
    match n {
        1 => x,
        0 => T::one(),
        _ if n >= 2 => powi(x, n%2) * powi(x*x, n/2),
        _ => T::one() / powi(x, -n)  // Here, we must have n < 0
    }
}
