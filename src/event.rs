//! This module defines the properties and storage of generated events

use crate::{
    linalg::{
        dimension::*,
        momentum::{Momentum, E},
        vecmat::*,
    },
    numeric::Float,
};

/// Number of incoming particles
pub const NUM_INCOMING: usize = 2;

/// Number of outgoing particles (replaces original INP)
pub const NUM_OUTGOING: usize = 3;

/// Number of particles in an event
pub const NUM_PARTICLES: usize = NUM_INCOMING + NUM_OUTGOING;

/// Vector type whose side is the number of particles in an event
pub type ParticleVector<T> = Vector5<T>;

/// Square matrix type whose side is the number of particles in an event
pub type ParticleMatrix<T> = Matrix5<T>;

/// Event data matrix type (columns are 4-coordinates, rows are particles)
type EventMatrix = Matrix5x4<Float>;

/// Slice of the event data matrix containing only outgoing particles
type OutgoingMomentaSlice<'matrix> = MatrixSlice<'matrix, Float, U3, U4, U1, U5>;

/// Row of the incoming electron in the event data matrix
pub const INCOMING_E_M: usize = 0;

/// Index of the incoming positron in the 4-momentum array
pub const INCOMING_E_P: usize = 1;

/// Number of possible spin values of the outgoing particles
pub const NUM_SPINS: usize = 2;

/// Storage for ee -> ppp event data
///
/// Encapsulates a vector of incoming and outgoing 4-momenta.
///
pub struct Event(EventMatrix);
//
impl Event {
    // ### CONSTRUCTOR ###

    /// Build an event from an event data matrix
    pub fn new(matrix: EventMatrix) -> Self {
        Self(matrix)
    }

    // ### ACCESSORS ###

    // FIXME: Remove a big bunch of those and use direct matrix access instead
    //        once const generics allow us to have elegant matrix slicing.

    /// Access the full internal 4-momentum array by reference
    pub fn all_momenta(&self) -> &EventMatrix {
        &self.0
    }

    /// Extract the 4-momentum of a single particle (internal for now)
    fn momentum(&self, par: usize) -> Momentum {
        Momentum::from_iterator(self.0.fixed_rows::<U1>(par).iter().cloned())
    }

    /// Extract the electron 4-momentum
    pub fn electron_momentum(&self) -> Momentum {
        self.momentum(INCOMING_E_M)
    }

    /// Extract the positron 4-momentum
    #[allow(dead_code)]
    pub fn positron_momentum(&self) -> Momentum {
        self.momentum(INCOMING_E_P)
    }

    /// Extract a photon's 4-momentum
    pub fn outgoing_momentum(&self, par: usize) -> Momentum {
        self.momentum(NUM_INCOMING + par)
    }

    /// Access the outgoing 4-momenta
    pub fn outgoing_momenta(&self) -> OutgoingMomentaSlice {
        self.0.fixed_rows::<U3>(NUM_INCOMING)
    }

    /// Minimal outgoing photon energy
    pub fn min_photon_energy(&self) -> Float {
        if cfg!(feature = "no-photon-sorting") {
            let first_out_e = self.outgoing_momenta()[(0, E)];
            self.outgoing_momenta()
                .fixed_columns::<U1>(E)
                .iter()
                .skip(1)
                .fold(first_out_e, |e1, &e2| if e1 < e2 { e1 } else { e2 })
        } else {
            // Use the fact that photons are sorted by decreasing energy
            self.outgoing_momenta()[(NUM_OUTGOING - 1, E)]
        }
    }

    // ### DEBUGGING TOOLS ###

    /// Dump 4-momenta of the 3 outgoing photons
    #[allow(dead_code)]
    pub fn display(&self) {
        let p_out = self.outgoing_momenta();
        for coord in 0..4 {
            print!("{}\t", coord);
            for part in 0..NUM_OUTGOING {
                print!("{}\t", p_out[(part, coord)]);
            }
            println!();
        }
    }
}
