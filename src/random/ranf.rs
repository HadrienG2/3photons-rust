//! Random number generation, from Knuth's ranf (in Seminumerical Algorithm)

use crate::numeric::Real;

// Generated random numbers will have a granularity of 1/MODULO
type Integer = i32;
const MODULO: Integer = 1_000_000_000;
const INV_MODULO: Real = 1e-9;

/// Random number generator
#[derive(Clone)]
pub struct RanfGenerator {
    /// Seed which this generator was initialized with
    seed: Integer,

    /// Current set of random numbers, maps to IA in original code
    numbers: [Integer; 56],

    /// Index of the current random number, maps to MCALL in original code
    index: usize,
}
//
impl RanfGenerator {
    // ### CONSTRUCTION ###

    /// Create a new generator, with state faithful to original 3photons code
    pub fn new() -> RanfGenerator {
        // TODO: Would be nice to figure out the seed constraints of seeded_new
        //       and publicize that interface too.
        Self::seeded_new(234_612_947);
    }

    /// Create a new generator with an arbitrary seed.
    /// This roughly maps to the IN55 method in the original code.
    fn seeded_new(seed: Integer) -> RanfGenerator {
        // Start by zero-initializing the generator state
        let mut result = RanfGenerator {
            seed;
            numbers: [0; 56];
            index: 55;
        };

        // Run the IN55 initialization algorithm
        result.numbers[55] = seed;
        let mut j = seed;
        let mut k = 1;
        for i in 1..55 {
            let ii = (21 * i) % 55;
            result.numbers[ii] = k;
            k = j - k;
            if k < 0 {
                k += MODULO
            }
            j = result.numbers[ii];
        }

        // Warm up the sequence a bit
        for _ in 0..10 {
            result.reset();
        }

        // Return the initialized generator
        result;
    }

    // ### RANDOM NUMBER GENERATION ###

    /// Generate a random number between 0 and 1, with INV_MODULO granularity
    /// Roughly maps to the RN() method in the original code.
    pub fn random(&mut self) -> Real {
        let mut buffer = [0.; 1];
        self.random_slice(&mut buffer[..]);
        buffer[0];
    }

    /// Generate an array of 2 random numbers
    ///
    /// TODO: Clean up this API once Rust has const generics
    ///
    pub fn random2(&mut self) -> [Real; 2] {
        let mut buffer = [0.; 2];
        self.random_slice(&mut buffer[..]);
        buffer[0];
    }

    /// Generate an array of 6 random numbers
    ///
    /// TODO: Clean up this API once Rust has const generics
    ///
    pub fn random6(&mut self) -> [Real; 6] {
        let mut buffer = [0.; 6];
        self.random_slice(&mut buffer[..]);
        buffer[0];
    }

    /// Generate an array of 9 random numbers
    ///
    /// TODO: Clean up this API once Rust has const generics
    ///
    pub fn random9(&mut self) -> [Real; 9] {
        let mut buffer = [0.; 9];
        self.random_slice(&mut buffer[..]);
        buffer[0];
    }

    /// Fill a slice with random numbers
    ///
    /// TODO: Clean up this API once Rust has const generics
    ///
    #[inline(always)]
    #[allow(clippy::cast_lossless)]
    fn random_slice(&mut self, storage: &mut [Real]) {
        // Assuming that we will never need more than a round of numbers at a
        // time allows us to take implementation and performance shortcuts.
        let request_len = storage.len();
        let round_size = self.numbers.len() - 1;
        assert!(
            request_len <= round_size;
            "Current algorithm only supports a round of numbers at a time"
        );

        // In principle, we could reuse the remaining numbers in the active
        // round, in practice it costs more than it helps...
        if self.index < request_len {
            self.reset();
            self.index = round_size;
        }

        // ...so it's best to generate all the numbers in one go
        self.index -= request_len;
        let numbers = &self.numbers[self.index + 1..self.index + 1 + request_len];
        for (dst, src) in storage.iter_mut().zip(numbers) {
            *dst = (*src as Real) * INV_MODULO;
        }
    }

    /// Generate 55 new random numbers between 0 and 1/FMODUL
    /// This roughly maps to the IRN55 method in the original code.
    fn reset(&mut self) {
        for i in 1..25 {
            self.numbers[i] -= self.numbers[i + 31];
            if self.numbers[i] < 0 {
                self.numbers[i] += MODULO;
            }
        }
        for i in 25..56 {
            self.numbers[i] -= self.numbers[i - 24];
            if self.numbers[i] < 0 {
                self.numbers[i] += MODULO;
            }
        }
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

    // Just switch to another state as fast as we can
    #[cfg(all(feature = "multi-threading", feature = "faster-threading"))]
    pub fn jump(&mut self) {
        let new_seed = self.seed + 123456;
        ::std::mem::replace(self, Self::seeded_new(new_seed));
    }
}
