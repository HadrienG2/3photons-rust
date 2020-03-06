//! This module defines the properties and storage of generated events

use crate::{
    linalg::{
        dimension::*,
        momentum::{Momentum, E},
        vecmat::*,
    },
    numeric::Float,
};
use std::fmt::Display;

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
    /// Build an event from an event data matrix
    pub fn new(matrix: EventMatrix) -> Self {
        Self(matrix)
    }

    // FIXME: Remove a big bunch of those and use direct matrix access instead
    //        once const generics allow us to have elegant matrix slicing.

    /// Access the full internal 4-momentum array by reference
    pub fn all_momenta(&self) -> &EventMatrix {
        &self.0
    }

    /// Extract the 4-momentum of a single particle (internal for now)
    ///
    /// We return an owned momentum vector because that's more convenient to
    /// handle than a slice and the compiler is smart enough to elide the copy.
    ///
    fn momentum(&self, par: usize) -> Momentum {
        Momentum::from_iterator(self.0.row(par).iter().copied())
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
            // This somewhat elaborate method (compared to, say, min_by())
            // allows us to propagate what we know about the matrix' layout and
            // contents to the compiler, enabling better optimization.
            let first_out_e = self.outgoing_momentum(0)[E];
            self.outgoing_momenta()
                .column(E)
                .iter()
                .skip(1)
                .fold(first_out_e, |e1, &e2| if e1 < e2 { e1 } else { e2 })
        } else {
            // Use the fact that photons are sorted by decreasing energy
            self.outgoing_momentum(NUM_OUTGOING - 1)[E]
        }
    }
}

impl Display for Event {
    /// Dump 4-momenta of the 3 outgoing photons
    fn fmt(&self, fmt: &mut std::fmt::Formatter<'_>) -> std::fmt::Result {
        let p_out = self.outgoing_momenta();
        for coord in 0..4 {
            write!(fmt, "{}\t", coord)?;
            for part in 0..NUM_OUTGOING {
                write!(fmt, "{}\t", p_out[(part, coord)])?;
            }
            writeln!(fmt)?;
        }
        Ok(())
    }
}
