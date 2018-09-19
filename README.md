# Rust version of 3photons

To run this version, make sure that you have a reasonably recent version of the
rustc compiler and the Cargo package manager around. You will find installation
instructions at <https://www.rust-lang.org/install.html>.

Once you're ready, the program can be built and run with the following command:

    cargo run --release


## Tuning knobs

By default, the simulation aims for maximal result reproducibility with respect
to the original `3photons` program. However, if you feel like going beyond that,
here are some features which you can enable:

- The `f32` feature moves all computations to single precision.
- The `fast-sincos` feature uses a different algorithm for computing points on
  the unit circle, which provides a great computational speedup at the cost of
  making reproducible parallelization more challenging.
- The `multi-threading` feature parallelizes the computation using multi-
  threading. By default, it generates results which are identical to the
  sequential version, but this has a performance, memory footprint and
  scalability cost. You can remove this constraint and allow the program to take
  more performance shortcuts by also enabling the `faster-threading` feature.
- The `no-photon-sorting` feature disables the sorting of outgoing photons by
  energy that used to be present in the original `3photons` code.
- The `standard-random` uses standard Rust abstractions for random number
  generation. The current algorithm (`xoshiro256+`) is not faster than what
  3photons uses by default (`ranf`), but is more friendly to parallelization
  (especially in `faster-threading` mode).

These features are enabled using the `--features` flag to cargo, as follows:

    cargo run --release --features multi-threading,faster-threading


## Reproducibility considerations

This version aims to produce results which are as close as possible to what the
original `3photons` program would emit, with one notable exception.

Instead of aggregating data about all events in a single accumulator, then
normalizing, it slices up the generated events into batches of a certain size
(implementation-defined, currently 10k events), accumulates them, and then
merges the resulting accumulators.

There are two good reasons to do this:

- It reduces worst-case accumulation error. As an accumulator grows large, there
  is a risk of new input not being well integrated into the accumulator.
- It makes the simulation more naturally parallelizable, enabling perfect
  reproducibility between the sequential and parallel versions.


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

As of 2018-09-17, various other things can be said about this version of
3photons, which I do not expect to remain accurate in the future.

1. You will need at least Rust v1.27 to compile it.
2. Comparing its performance with that of other versions would be unfair,
   because the porting process enabled many optimizations which other versions
   did not receive yet.
