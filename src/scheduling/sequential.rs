//! Sequential back-end of the simulation

use crate::{random::RandomGenerator, resfin::ResultsAccumulator, scheduling::EVENT_BATCH_SIZE};

/// Simulate events in sequential mode
///
/// We use batched logic even in sequential mode, in order to achieve
/// reproducibility with respect to multi-threaded runs.
///
/// Note that this is anyways generally a good thing to do when accumulating
/// lots of results, as otherwise the accumulator will eventually grow much
/// larger than the accumulated values and numerical accumulation errors
/// will start to blow up.
///
pub fn run_simulation_impl<'cfg>(
    mut num_events: usize,
    mut rng: RandomGenerator,
    simulate_events: impl Send + Sync + Fn(usize, &mut RandomGenerator) -> ResultsAccumulator<'cfg>,
) -> ResultsAccumulator<'cfg> {
    // Some double-checking cannot hurt...
    assert!(num_events > 0, "Must simulate at least one event");

    // Initialize the accumulator with the first batch of events
    let first_batch_size = ::std::cmp::min(EVENT_BATCH_SIZE, num_events);
    num_events -= first_batch_size;
    let mut accumulator = simulate_events(first_batch_size, &mut rng);

    // Simulate and integrate complete batches of events (if any)
    let num_full_batches = num_events / EVENT_BATCH_SIZE;
    for _ in 0..num_full_batches {
        accumulator.merge(simulate_events(EVENT_BATCH_SIZE, &mut rng));
    }
    num_events %= EVENT_BATCH_SIZE;

    // Integrate the remaining events
    accumulator.merge(simulate_events(num_events, &mut rng));

    // Return the final accumulated results
    accumulator
}
