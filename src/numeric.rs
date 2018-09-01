//! Basic numerical concepts used throughout the program

use num_complex;


/// Floating-point precision is configured here
pub use std::f64 as reals;
pub type Real = f64;
pub type Complex = num_complex::Complex<Real>;

/// Mathematical functions
pub mod functions {
    use num_traits::One;
    use std::ops::{Div, Mul};
    use super::Real;


    // Generate trivial wrappers around member functions of floating-point types
    // to allow for a consistent notation in floating-point expressions. For
    // example, this allows using sin(x) instead of x.sin().
    //
    macro_rules! prefix_unary_funcs {
        ( $( $name:ident ),* ) => ( $(
            #[inline]
            pub fn $name(x: Real) -> Real {
                x.$name()
            }
        )* )
    }
    //
    prefix_unary_funcs!{ abs, cos, exp, ln, sin, sqrt }


    // The following functions are generic in nature, and designed to work on
    // both real and complex numbers.

    /// Compute the square of any number, optimized shorthand for powi(x, 2)
    #[inline]
    pub fn sqr<T>(x: T) -> T
        where T: Mul<Output=T> + Copy
    {
        x * x
    }

    /// Raise a number to an arbitrary integer power, like {float}.powi(), using
    /// a binary exponentiation algorithm.
    ///
    /// This recursive version is clean and concise, but should only be used
    /// when n is known at compile time, and the recursion can be unrolled...
    ///
    #[inline]
    pub fn powi<T>(x: T, n: i32) -> T
        where T: Mul<Output=T> + Div<Output=T> + Copy + One
    {
        match n {
            0 => T::one(),
            1 => x,
            _ if n >= 2 => powi(x, n%2) * powi(sqr(x), n/2),
            _ => T::one() / powi(x, -n)  // Here, we must have n < 0
        }
    }
}