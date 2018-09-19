//! Random number generation module. Uses either a port of 3photon's ranf random
//! number generator or the "rand" crate that is the Rust standard for RNGs.

/// Build the ranf module if it's going to be used
#[cfg(not(feature = "standard-random"))]
mod ranf;

/// Else mimick the interface of RanfGenerator on top of the rand crate
#[cfg(feature = "standard-random")]
mod standard_random {
    use ::numeric::Real;
    use rand::Rng;

    /// Random number generation engine in use
    type Engine = ::xoshiro::Xoshiro256Plus;

    /// Facade which makes the rand crate look like RanfGenerator
    #[derive(Clone)]
    pub struct RandGenerator {
        rng: Engine,
    }
    //
    impl RandGenerator {
        // Spawn a new random number generator
        pub fn new() -> Self {
            Self {
                rng: Engine::from_seed_u64(12345),
            }
        }

        // Generate a random floating-point number between 0 and 1
        pub fn random(&mut self) -> Real {
            self.rng.gen()
        }

        /// Generate an array of 2 random numbers
        ///
        /// TODO: Clean up this API once Rust has const generics
        ///
        pub fn random2(&mut self) -> [Real; 2] {
            self.rng.gen()
        }

        /// Generate an array of 3 random numbers
        ///
        /// TODO: Clean up this API once Rust has const generics
        ///
        pub fn random3(&mut self) -> [Real; 3] {
            self.rng.gen()
        }

        /// Generate an array of 6 random numbers
        ///
        /// TODO: Clean up this API once Rust has const generics
        ///
        pub fn random6(&mut self) -> [Real; 6] {
            self.rng.gen()
        }

        /// Generate an array of 9 random numbers
        ///
        /// TODO: Clean up this API once Rust has const generics
        ///
        pub fn random9(&mut self) -> [Real; 9] {
            self.rng.gen()
        }

        // Advance state as if "iterations" numbers had been generated
        #[cfg(all(feature = "multi-threading",
                  not(feature = "faster-threading")))]
        pub fn skip(&mut self, iterations: usize) {
            for _ in 0..iterations { self.random(); }
        }

        // Advance state in an arbitrary but maximally fast way
        #[cfg(all(feature = "multi-threading",
                  feature = "faster-threading"))]
        pub fn jump(&mut self) {
            self.rng.jump();
        }
    }
}

/// Select the RNG implementation in use
#[cfg(feature = "standard-random")]
pub use self::standard_random::RandGenerator as RandomGenerator;
#[cfg(not(feature = "standard-random"))]
pub use self::ranf::RanfGenerator as RandomGenerator;