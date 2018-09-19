//! Random number generation, from Knuth's ranf (in Seminumerical Algorithm)

use ::numeric::Real;

use std::cmp;


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
        // FIXME: Provide control over seeding
        Self::seeded_new(234612947)
    }

    /// Create a new generator with an arbitrary seed.
    /// This roughly maps to the IN55 method in the original code.
    fn seeded_new(seed: Integer) -> RanfGenerator {
        // Start by zero-initializing the generator state
        let mut result = RanfGenerator {
            seed,
            numbers: [0; 56],
            index: 55,
        };

        // Run the IN55 initialization algorithm
        result.numbers[55] = seed;
        let mut j = seed;
        let mut k = 1;
        for i in 1..55 {
            let ii = (21 * i) % 55;
            result.numbers[ii] = k;
            k = j - k;
            if k < 0 { k += MODULO };
            j = result.numbers[ii];
        }

        // Warm up the sequence a bit
        for _ in 0..10 {
            result.reset();
        }

        // Return the initialized generator
        result
    }


    // ### RANDOM NUMBER GENERATION ###

    /// Generate a random number between 0 and 1, with INV_MODULO granularity
    /// Roughly maps to the RN() method in the original code.
    pub fn random(&mut self) -> Real {
        let mut buffer = [0.; 1];
        self.random_slice(&mut buffer[..]);
        buffer[0]
    }

    /// Generate an array of 2 random numbers
    ///
    /// TODO: Clean up this API once Rust has const generics
    ///
    pub fn random2(&mut self) -> [Real; 2] {
        let mut buffer = [0.; 2];
        self.random_slice(&mut buffer[..]);
        buffer
    }

    /// Generate an array of 9 random numbers
    ///
    /// TODO: Clean up this API once Rust has const generics
    ///
    pub fn random9(&mut self) -> [Real; 9] {
        let mut buffer = [0.; 9];
        self.random_slice(&mut buffer[..]);
        buffer
    }

    /// Fill a slice with random numbers
    ///
    /// TODO: Clean up this API once Rust has const generics
    ///
    #[inline(always)]
    fn random_slice(&mut self, storage: &mut [Real]) {
        // Implementation simplification
        assert!(storage.len() <= 54,
                "Current algorithm only supports 54 numbers at a time");

        // Generate as much as possible without an RNG reset
        let first_chunk = cmp::min(self.index, storage.len());
        let (storage1, storage2) = storage.split_at_mut(first_chunk);
        self.index -= storage1.len();
        storage1.iter_mut().zip(&self.numbers[self.index+1..])
                .for_each(|(dst, src)| *dst = (*src as Real) * INV_MODULO);

        // Reset the RNG (if needed) and generate remaining numbers (if any)
        if self.index == 0 {
            self.reset();
            self.index = 55 - storage2.len();
            storage2.iter_mut().zip(&self.numbers[self.index+1..])
                    .for_each(|(dst, src)| *dst = (*src as Real) * INV_MODULO);
        }
    }


    /// Generate 55 new random numbers between 0 and 1/FMODUL
    /// This roughly maps to the IRN55 method in the original code.
    fn reset(&mut self) {
        for i in 1..25 {
            self.numbers[i] -= self.numbers[i+31];
            if self.numbers[i] < 0 { self.numbers[i] += MODULO };
        }
        for i in 25..56 {
            self.numbers[i] -= self.numbers[i-24];
            if self.numbers[i] < 0 { self.numbers[i] += MODULO };
        }
    }

    // Advance state as if "iterations" numbers had been generated
    #[cfg(all(feature = "multi-threading",
              not(feature = "faster-threading")))]
    pub fn skip(&mut self, iterations: usize) {
        for _ in 0..iterations { self.random(); }
    }

    // Just switch to another state as fast as we can
    #[cfg(all(feature = "multi-threading",
              feature = "faster-threading"))]
    pub fn jump(&mut self) {
        let new_seed = self.seed + 123456789;
        ::std::mem::replace(self, Self::seeded_new(new_seed));
    }
}
