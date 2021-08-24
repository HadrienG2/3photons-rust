//! Random number generation, from Knuth's ranf (in Seminumerical Algorithm)

use crate::numeric::Float;

// Generated random numbers will have a granularity of 1/MODULO
type RanfInt = i32;
const MODULO: RanfInt = 1_000_000_000;
const INV_MODULO: Float = 1e-9;

/// Random number generator
#[derive(Clone)]
pub struct RanfGenerator {
    /// Seed which this generator was initialized with
    seed: RanfInt,

    /// Current set of random numbers, maps to IA in original code
    numbers: [RanfInt; 56],

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
        Self::seeded_new(234_612_947)
    }

    /// Create a new generator with an arbitrary seed.
    /// This roughly maps to the IN55 method in the original code.
    fn seeded_new(seed: RanfInt) -> RanfGenerator {
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
        result
    }

    // ### RANDOM NUMBER GENERATION ###

    /// Generate a random number between 0 and 1, with INV_MODULO granularity
    /// Roughly maps to the RN() method in the original code.
    pub fn random(&mut self) -> Float {
        self.random_array::<1>()[0]
    }

    /// Generate an array of random numbers
    #[allow(clippy::cast_lossless)]
    pub fn random_array<const N: usize>(&mut self) -> [Float; N] {
        // Assuming that we will never need more than a round of numbers at a
        // time allows us to take implementation and performance shortcuts.
        let round_size = self.numbers.len() - 1;
        assert!(
            N <= round_size,
            "Current algorithm only supports a round of numbers at a time"
        );

        // In principle, we could reuse the remaining numbers in the active
        // round, in practice it costs more than it helps...
        if self.index < N {
            self.reset();
            self.index = round_size;
        }

        // ...so it's best to generate all the numbers in one go
        self.index -= N;
        let numbers = &self.numbers[self.index + 1..self.index + 1 + N];
        let mut result = [0.; N];
        for (dst, src) in result.iter_mut().zip(numbers) {
            *dst = (*src as Float) * INV_MODULO;
        }
        result
    }

    /// Generate 55 new random numbers between 0 and 1/FMODUL
    /// This roughly maps to the IRN55 method in the original code.
    fn reset(&mut self) {
        for i in 1..25 {
            self.numbers[i] -= self.numbers[i + 31];
            if self.numbers[i] < 0 {
                self.numbers[i] += MODULO
            }
        }
        for i in 25..56 {
            self.numbers[i] -= self.numbers[i - 24];
            if self.numbers[i] < 0 {
                self.numbers[i] += MODULO
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

    // Advance state as if random_array() had been called
    #[cfg(all(feature = "multi-threading", not(feature = "faster-threading")))]
    pub fn skip_array<const N: usize>(&mut self) {
        self.random_array::<N>();
    }

    // Just switch to another state as fast as we can
    #[cfg(all(feature = "multi-threading", feature = "faster-threading"))]
    pub fn jump(&mut self) {
        let new_seed = self.seed + 123456;
        *self = Self::seeded_new(new_seed);
    }
}
