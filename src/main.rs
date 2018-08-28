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

// `error_chain!` can recurse deeply
#![recursion_limit = "1024"]

#[macro_use] extern crate error_chain;

extern crate chrono;
extern crate nalgebra;
extern crate num_complex;
extern crate num_traits;

mod config;
mod coupling;
mod evcut;
mod event;
mod linalg;
mod numeric;
mod output;
mod ranf;
mod rescont;
mod resfin;
mod scalar;
mod spinor;

use config::Configuration;
use coupling::Couplings;
use event::Event;
use ranf::RanfGenerator;
use rescont::ResultContribution;
use resfin::ResultsBuilder;
use scalar::ScalarProducts;
use spinor::SpinorProducts;
use std::time::Instant;


/// This will act as our main function, with suitable error handling
quick_main!(|| {
    // ### CONFIGURATION READOUT ###

    // The work of loading, parsing, and checking the configuration has now been
    // offloaded to a dedicated struct
    let cfg = Configuration::new("valeurs")
                            .chain_err(|| "Failed to load the configuration")?;


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

    // Initialize the random number generator
    let mut rng = RanfGenerator::new();

    // Initialize the event generator
    let mut event = Event::new(cfg.e_tot);

    // Initialize results accumulator
    let mut res_builder = ResultsBuilder::new(&cfg, event.weight());


    // ### SAMPLING LOOP ###

    for _ in 0..cfg.num_events {
        // Generate an event
        event.generate_momenta(&mut rng);

        // Sort the outgoing photons by energy
        event.sort_output_momenta();

        // Compute spinor inner products and scalar products
        let spinor = SpinorProducts::new(&event);
        let scalar = ScalarProducts::new(&spinor);

        // If the event passes the cut, compute the total weight (incl. matrix
        // elements) and integrate it into the final results.
        if cfg.event_cut.keep(&event, scalar) {
            let res_contrib = ResultContribution::new(&couplings, &spinor);
            // NOTE: This is where the original code would display the result
            res_builder.integrate(res_contrib);
            // NOTE: This is where the FORTRAN code would fill histograms
        }
    }


    // ### PHYSICAL RESULTS COMPUTATION ###
    
    // Now that we have integrated everything, compute the physical results
    let res_fin = res_builder.finalize();

    // NOTE: This is where the FORTRAN code would normalize histograms


    // ### RESULTS DISPLAY AND STORAGE ###

    // Measure how much time has elapsed
    let elapsed_time = saved_time.elapsed();
    
    // Send the results to the standard output and to disk and we're done
    output::dump_results(&cfg, res_fin, elapsed_time)
           .chain_err(|| "Failed to output the results")
});



// Here are the various things that can go wrong during the main function
mod errors {
    error_chain!{
        links{
            // Something bad happened while loading the configuration
            Config(::config::Error, ::config::ErrorKind);
        }

        foreign_links{
            // Something bad happened while outputting the results
            ResultsIo(::std::io::Error);
        }
    }
}
//
use errors::*;
