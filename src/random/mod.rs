//! Random number generation module. Uses either a port of 3photon's ranf random
//! number generator or the "rand" crate that is the Rust standard for RNGs.

/// Build the ranf module if it's going to be used
#[cfg(not(feature = "standard-random"))]
mod ranf;

/// Else mimick the interface of RanfGenerator on top of the rand crate
#[cfg(feature = "standard-random")]
mod standard;

// Select the RNG implementation in use
#[cfg(not(feature = "standard-random"))]
pub use self::ranf::RanfGenerator as RandomGenerator;
#[cfg(feature = "standard-random")]
pub use self::standard::RandGenerator as RandomGenerator;
