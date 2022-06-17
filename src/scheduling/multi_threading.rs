//! Multi-threaded back-end of the simulation

#[cfg(feature = "faster-threading")]
use std::sync::atomic::{AtomicBool, Ordering};

use crate::{random::RandomGenerator, resacc::ResultsAccumulator, scheduling::EVENT_BATCH_SIZE};

use std::sync::Mutex;

/// Simulate events in multi-threaded mode
///
/// Takes the random number generator state as input, and does not finalize
/// the output results, so should be readily amenable to extra layers of
/// parallelization (such as distribution across multiple compute nodes).
///
pub fn run_simulation_impl<'cfg>(
    mut num_events: usize,
    mut rng: RandomGenerator,
    simulate_events: impl Send + Sync + Fn(usize, &mut RandomGenerator) -> ResultsAccumulator<'cfg>,
) -> ResultsAccumulator<'cfg> {
    // Some double-checking cannot hurt...
    assert!(num_events > 0, "Must simulate at least one event");

    // We know in advance how many batches of event we will process
    let num_batches = num_events / EVENT_BATCH_SIZE
        + if num_events % EVENT_BATCH_SIZE == 0 {
            0
        } else {
            1
        };

    // The results of parallel tasks will be aggregated...
    let accumulator = {
        // ...in a way that is optimized for numerical reproduciblity
        #[cfg(not(feature = "faster-threading"))]
        {
            ReproducibleAccumulator::new(num_batches)
        }

        // ...in a way that is optimized for computational performance
        #[cfg(feature = "faster-threading")]
        {
            FastAccumulator::new(num_batches)
        }
    };

    // This function is a synchronization scope: it will only return
    // once all inner tasks have been executed
    rayon::scope(|scope| {
        // For each requested batch of events...
        for batch_id in 0..num_batches {
            let batch_size = ::std::cmp::min(num_events, EVENT_BATCH_SIZE);
            num_events -= batch_size;

            // Spawn a task which is responsible for simulating them
            let mut task_rng = rng.clone();
            let accumulator_ref = &accumulator;
            let simulate_events_ref = &simulate_events;
            scope.spawn(move |_| {
                let result = simulate_events_ref(batch_size, &mut task_rng);
                accumulator_ref.set_task_result(batch_id, result);
            });

            // By default, we rigorously simulate the impact of each
            // task on the random number generator. This makes the
            // results reproducible, but slows down the scheduling
            // thread and may thus reduce performance and scalability.
            #[cfg(not(feature = "faster-threading"))]
            crate::evgen::EventGenerator::simulate_event_batch(&mut rng, batch_size);

            // In non-reproducible mode, we instead ask the RNG to
            // switch to a wildly different state as quickly as it can.
            #[cfg(feature = "faster-threading")]
            rng.jump();
        }
    });

    // Extract the results from the accumulator
    accumulator.get_merged_result()
}

/// Reproducibility-optimized results accumulation mechanism
#[cfg(not(feature = "faster-threading"))]
struct ReproducibleAccumulator<'cfg> {
    /// Storage for the intermediary simulation results of parallel tasks
    results: Box<[Mutex<Option<ResultsAccumulator<'cfg>>>]>,
}
//
#[cfg(not(feature = "faster-threading"))]
impl<'cfg> ReproducibleAccumulator<'cfg> {
    /// Set up results storage for N parallel tasks
    fn new(num_tasks: usize) -> Self {
        assert!(num_tasks > 0, "There should be at least one task");
        Self {
            results: (0..num_tasks)
                .map(|_| Mutex::new(None))
                .collect::<Vec<_>>()
                .into_boxed_slice(),
        }
    }

    /// Integrate the results of the n-th simulation task
    fn set_task_result(&self, task_id: usize, result: ResultsAccumulator<'cfg>) {
        let mut lock = self.results[task_id]
            .lock()
            .expect("Mutex data should be valid");
        assert!(lock.is_none(), "Tasks should not report results twice");
        *lock = Some(result);
    }

    /// Aggregate the results in a reproducible fashion
    fn get_merged_result(self) -> ResultsAccumulator<'cfg> {
        // Start iterating over the task results
        let mut results_iter = self.results.into_vec().into_iter().map(|entry| {
            entry
                .into_inner()
                .expect("Mutex data should be valid")
                .expect("Result should be ready")
        });

        // Initialize results storage with the result of the first task
        let first_result = results_iter
            .next()
            .expect("There should be at least one task");

        // Merge the results of the other tasks
        results_iter.fold(first_result, |mut r1, r2| {
            r1.merge(r2);
            r1
        })
    }
}

/// Speed-optimized results accumulation mechanism
#[cfg(feature = "faster-threading")]
struct FastAccumulator<'cfg> {
    /// Storage location in which results will be merged out of order
    merged_result: Mutex<Option<ResultsAccumulator<'cfg>>>,

    /// Truth that each task has reported its results
    task_finished: Box<[AtomicBool]>,
}
//
#[cfg(feature = "faster-threading")]
impl<'cfg> FastAccumulator<'cfg> {
    /// Set up results storage for N parallel tasks
    fn new(num_tasks: usize) -> Self {
        assert!(num_tasks > 0, "There should be at least one task");
        Self {
            merged_result: Mutex::new(None),
            task_finished: (0..num_tasks)
                .map(|_| AtomicBool::new(false))
                .collect::<Vec<_>>()
                .into_boxed_slice(),
        }
    }

    /// Integrate the results of the n-th simulation task
    #[allow(unknown_lints, clippy::significant_drop_in_scrutinee)]
    fn set_task_result(&self, task_id: usize, result: ResultsAccumulator<'cfg>) {
        // Initialize the accumulator or merge the task result into it
        match *self
            .merged_result
            .lock()
            .expect("Mutex data should be valid")
        {
            // If we are the first, initialize the accumulator
            ref mut storage @ None => *storage = Some(result),

            // Otherwise, merge our results with those that are already here
            Some(ref mut accumulator) => accumulator.merge(result),
        }

        // Remember that this task has completed its work
        let was_finished = self.task_finished[task_id].swap(true, Ordering::Relaxed);
        assert!(!was_finished, "Tasks should not set their result twice");
    }

    /// Aggregate the results in a reproducible fashion
    fn get_merged_result(self) -> ResultsAccumulator<'cfg> {
        // Check that all tasks have completed their work
        for ready in self.task_finished.into_vec().into_iter() {
            assert!(
                ready.load(Ordering::Relaxed),
                "All tasks should have completed their work"
            );
        }

        // Collect the merged result
        self.merged_result
            .into_inner()
            .expect("Mutex data should be valid")
            .expect("Result should be ready")
    }
}
