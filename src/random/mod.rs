//! Random number generation module. Uses either a port of 3photon's ranf random
//! number generator or the "rand" crate that is the Rust standard for RNGs.

#[cfg(not(feature = "standard-random"))]
mod ranf;

/// Mimick the interface of RanfGenerator on top of the rand crate
#[cfg(feature = "standard-random")]
mod standard_random {
    use ::numeric::Real;
    use rand::Rng;

    /// Random number generation engine in use
    type Engine = ::xoshiro::Xoshiro256Plus;

    /// Facade which makes the rand crate look like RanfGenerator
    pub struct RandGenerator {
        rng: Engine,
    }
    //
    impl RandGenerator {
        pub fn new() -> Self {
            Self {
                rng: Engine::from_seed_u64(12345),
            }
        }

        pub fn random(&mut self) -> Real {
            self.rng.gen()
        }
    }
}

/// Select the RNG implementation in use
#[cfg(feature = "standard-random")]
pub use self::standard_random::RandGenerator as RandomGenerator;
#[cfg(not(feature = "standard-random"))]
pub use self::ranf::RanfGenerator as RandomGenerator;