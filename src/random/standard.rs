//! This module mimicks the interface of RanfGenerator on top of the
//! abstractions of the standard "random" crate.

use crate::numeric::Float;
use rand::{Rng, SeedableRng};

// Select random number generation engine in use
#[cfg(feature = "f32")]
type Engine = rand_xoshiro::Xoshiro128Plus;
#[cfg(not(feature = "f32"))]
type Engine = rand_xoshiro::Xoshiro256Plus;

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
            rng: Engine::seed_from_u64(12345),
        }
    }

    // Generate a random floating-point number between 0 and 1
    pub fn random(&mut self) -> Float {
        self.rng.gen()
    }

    /// Generate an array of 2 random numbers
    ///
    /// TODO: Clean up this API once Rust has const generics
    ///
    pub fn random2(&mut self) -> [Float; 2] {
        self.rng.gen()
    }

    /// Generate an array of 6 random numbers
    ///
    /// TODO: Clean up this API once Rust has const generics
    ///
    pub fn random6(&mut self) -> [Float; 6] {
        self.rng.gen()
    }

    /// Generate an array of 9 random numbers
    ///
    /// TODO: Clean up this API once Rust has const generics
    ///
    pub fn random9(&mut self) -> [Float; 9] {
        self.rng.gen()
    }

    // Advance state as if random() had been called "iteration" times
    #[cfg(all(feature = "multi-threading", not(feature = "faster-threading")))]
    pub fn skip(&mut self, iterations: usize) {
        for _ in 0..iterations {
            self.random();
        }
    }

    // Advance state as if random9() had been called
    ///
    /// TODO: Clean up this API once Rust has const generics
    ///
    #[cfg(all(feature = "multi-threading", not(feature = "faster-threading")))]
    pub fn skip9(&mut self) {
        self.random9();
    }

    // Advance state in an arbitrary but maximally fast way
    #[cfg(all(feature = "multi-threading", feature = "faster-threading"))]
    pub fn jump(&mut self) {
        self.rng.jump();
    }
}
