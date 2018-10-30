//! This module takes care of scheduling the simulation work, encapsulating use
//! of multiple threads and anything else that will come in the future

#[cfg(not(feature = "multi-threading"))] mod sequential;
#[cfg(feature = "multi-threading")] mod multi_threading;

use crate::{
    random::RandomGenerator,
    resfin::{FinalResults, ResultsBuilder},
};


/// Size of the simulated event batches
///
/// Simulated events are grouped in batches of a certain size in order to reduce
/// accumulation error and achieve perfect reproducibility between sequential
/// and parallel runs of the simulation.
///
/// This constant may need to be tuned in the future if CPUs become faster or
/// synchronization overhead changes. But the rate of such change is expected to
/// be low enough for hard-coding of this constant to be reasonable.
///
const EVENT_BATCH_SIZE: usize = 10_000;


/// Run the simulation in the manner that was configured at build time.
///
/// Takes as parameters the total number of events to be simulated, and a
/// simulation kernel that simulates a certain number of events given an initial
/// random number generator state.
///
/// Returns the finalized simulation results
///
pub fn run_simulation<'a>(
    num_events: usize,
    simulate_events: impl Send + Sync
                          + Fn(usize,
                               &mut RandomGenerator) -> ResultsBuilder<'a>
) -> FinalResults<'a> {
    // Check that the user is being reasonable (should have already been checked
    // at configuration time, but bugs can happen...)
    assert!(num_events > 0, "Must simulate at least one event");

    // Initialize the random number generator
    let rng = RandomGenerator::new();

    // Integrate simulation results...
    let results_builder = {
        // ...in sequential mode
        #[cfg(not(feature = "multi-threading"))]
        { sequential::run_simulation_impl(num_events,
                                          rng,
                                          simulate_events) }

        // ...in multi-threaded mode
        #[cfg(feature = "multi-threading")]
        { multi_threading::run_simulation_impl(num_events,
                                               rng,
                                               simulate_events) }
    };

    // Finalize the results
    results_builder.finalize()
}