//! 3 photons: A simple Monte Carlo simulation
//!
//!
//! # Introduction (for the physicist)
//!
//! This small computational program computes cross-section for the particle
//! physics process electron + positron gives three photons (e⁺e⁻ → 𝛾𝛾𝛾).
//!
//! It distinguishes a classical Standard Model contribution, of purely Quantum
//! ElectroDynamic origin and an hypothetic, beyond the Standard Model, New
//! Physics contribution, phenomenologically described by two effective
//! operators.
//!
//! It was designed in the LEP era, so these new interactions occurs between the
//! Z⁰ boson and the three photons.
//!
//! The effective operator can be related to specific models, among which
//! magnetic monopoles that run in a four points loop. The two operators exhibit
//! different (???)
//!
//!
//! # Introduction (for the numerical guy)
//!
//! The physicist want to compute a (multidimensional) integral, so we chose a
//! Monte Carlo algorithm
//!
//!
//! # Introduction (for the computer guy)
//!
//! this program started in a purely procedural style:
//!
//! * read in parameters and initialise counters
//! * loop over (random) event,
//!     * determining their geometrical and energy configuration,
//!     * their phase space weight,
//!     * their transition probability for each polarisation/helicity
//!       configuration, depending on coupling strength
//!     * sum it up
//! * then display / store the result.
//!
//! The use of common (for the original Fortran) or struct (in C) or record
//! types (in Ada) or classes (in C++) illustrates an object oriented design.
//!
//! The fact that we can plug each phase's output as the input of the next phase
//! lend to a functionnal approach.

#![warn(missing_docs)]

mod config;
mod coupling;
mod evcut;
mod event;
mod linalg;
mod numeric;
mod output;
mod random;
mod rescont;
mod resfin;
mod scheduling;
mod spinor;

use anyhow::Context;

use crate::{
    config::Configuration, coupling::Couplings, event::EventGenerator, random::RandomGenerator,
    rescont::ResultContribution, resfin::ResultsBuilder,
};

use std::time::Instant;

/// We'll use anyhow's type-erased result type throughout the application
type Result<T> = anyhow::Result<T>;

/// This will act as our main function, with suitable error handling
fn main() -> Result<()> {
    // ### CONFIGURATION READOUT ###

    // The work of loading, parsing, and checking the configuration has now been
    // offloaded to a dedicated struct
    let cfg = Configuration::load("valeurs").context("Failed to load the configuration")?;

    // ### SIMULATION INITIALIZATION ###

    // NOTE: Unlike the C++ version, we start the clock after configuration I/O,
    //       to avoid IO-induced timing fluctuations
    let saved_time = Instant::now();

    // NOTE: Removed final particle mass array. Since we are simulating photons,
    //       we know the masses to be zero at compile time.

    // NOTE: Deleted the original WTEV value. In the C++ code, it was
    //       overwritten by the first RAMBO call w/o having ever been read!

    // Compute physical couplings
    let couplings = Couplings::new(&cfg);

    // Initialize the event generator
    let evgen = EventGenerator::new(cfg.e_tot);

    // ### SIMULATION EXECUTION ###

    // This kernel simulates a number of events, given an initial random number
    // generator state, and return the accumulated intermediary results
    let simulate_events = |num_events: usize, rng: &mut RandomGenerator| -> ResultsBuilder {
        // Setup a results accumulator
        let mut res_builder = ResultsBuilder::new(&cfg, evgen.event_weight());

        // Simulate the requested number of events
        for _ in 0..num_events {
            // Generate an event
            let event = evgen.generate(rng);

            // If the event passes the cut, compute the total weight (incl.
            // matrix elements) and integrate it into the final results.
            if cfg.event_cut.keep(&event) {
                let res_contrib = ResultContribution::new(&couplings, &event);
                // NOTE: The original code would display the result here
                res_builder.integrate(res_contrib);
                // NOTE: The FORTRAN code would fill histograms here
            }
        }

        // Return the accumulated results
        res_builder
    };

    // Run the simulation
    let result = scheduling::run_simulation(cfg.num_events, simulate_events);

    // NOTE: This is where the FORTRAN code would normalize histograms

    // ### RESULTS DISPLAY AND STORAGE ###

    // Measure how much time has elapsed
    let elapsed_time = saved_time.elapsed();

    // Send the results to the standard output and to disk and we're done
    output::dump_results(&cfg, &result, elapsed_time).context("Failed to output the results")?;

    // ...and we're done
    Ok(())
}
