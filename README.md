# Rust version of 3photons

To run this version, make sure that you have a reasonably recent version of the
rustc compiler and the Cargo package manager around. You will find installation
instructions at <https://www.rust-lang.org/install.html>.

Once you're ready, the program can be built and run with the following command:

    cargo run --release


## Tuning knobs

By default, the simulation aims for ~perfect result reproducibility with respect
to the original 3photons program. However, if you feel like exploring other
points in the design space, here are some features which you can enable:

- The `f32` feature moves all computations to single precision.
- The `fast-sincos` feature uses a different algorithm for computing points on
  the unit circle, which provides a great computational speedup.
- The `standard-random` uses standard Rust abstractions for random number
  generation. The current algorithm (`xoshiro256+`) is _slower_ than what
  3photons uses by default (`ranf`), but more amenable to parallelization.

To enable all of these feautres, one would run:

    cargo run --release --features fast-sincos,f32,standard-random


## Error handling

Owing to Rust's strong focus on software engineering best practices, this
version received more thoughts than most others in various areas including
immutability by default, floating-point resolution independence, global variable
avoidance, and error handling.

On this later front, a key design goal is that all known error conditions should
abort the program, and most of them will output a clear error message when doing
so. Failure to meet this standard should be considered as a bug and reported.

Where a compromise between performance and error detection exists, the
associated error checking is only enabled in debug mode. To run the program
in this configuration, use the following command:

    cargo run


## Miscellaneous remarks, subjected to future obsolescence

As of 2017-04-20, various other things can be said about this version of
3photons, which I do not expect to remain accurate in the future.

1. You will need at least Rust v1.16 to compile it.
2. Comparing its performance with that of other versions would be unfair,
   because the porting process enabled many optimizations which other versions
   did not receive yet.
