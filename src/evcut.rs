//! Mechanism to apply a cut to generated events

use crate::{
    event::{Event, NUM_OUTGOING},
    momentum::{E, X},
    numeric::Float,
};
use prefix_num_ops::real::*;

/// Cuts on generated events
pub struct EventCut {
    /// Cut on maximum cosine of (beam, photons) angle
    pub beam_photons_cut: Float,

    /// Cut on maximum cosine of (photon, photon) angle
    pub photon_photon_cut: Float,

    /// Cut on minimum photon energy
    pub e_min: Float,

    /// Cut on minimum cosine of (beam, normal to the photon plane) angle
    pub beam_photon_plane_cut: Float,
}
//
impl EventCut {
    /// Setup the cuts on generated events
    pub fn new(
        beam_photons_cut: Float,
        photon_photon_cut: Float,
        e_min: Float,
        beam_photon_plane_cut: Float,
    ) -> Self {
        EventCut {
            beam_photons_cut,
            photon_photon_cut,
            e_min,
            beam_photon_plane_cut,
        }
    }

    /// Decide whether a generated event passes the cut or should be rejected
    pub fn keep(&self, event: &Event) -> bool {
        // Check if the outgoing photons pass the energy cut
        if event.min_photon_energy() < self.e_min {
            return false;
        }

        // Get the incoming electron 4-momentum
        let p_el = event.electron_momentum();

        // Check if the (beam, photon) angles pass the cut
        {
            let ps_out = event.outgoing_momenta();
            let ps_out_xyz = ps_out.fixed_columns::<3>(X);
            let cos_nums = ps_out_xyz * p_el.xyz();
            let cos_denoms = ps_out.column(E) * p_el[E];
            for (&num, denom) in cos_nums.iter().zip(cos_denoms.iter()) {
                if abs(num) > self.beam_photons_cut * denom {
                    return false;
                }
            }
        }

        // Check if the (photon1, photon{2, 3}) angles pass the cut
        for ph1 in 0..NUM_OUTGOING - 1 {
            // FIXME: Consider re-vectorizing this inner loop once const
            //        generics make it ergonomic to do so.
            for ph2 in ph1 + 1..NUM_OUTGOING {
                let p_ph1 = event.outgoing_momentum(ph1);
                let p_ph2 = event.outgoing_momentum(ph2);
                let cos_num = p_ph1.xyz().dot(&p_ph2.xyz());
                let cos_denom = p_ph1[E] * p_ph2[E];
                if cos_num > self.photon_photon_cut * cos_denom {
                    return false;
                }
            }
        }

        // Compute a vector which is normal to the outgoing photon plane
        // This notion is only valid because we have three output photons
        assert_eq!(NUM_OUTGOING, 3, "This part assumes 3 outgoing particles");
        let n_ppp = event
            .outgoing_momentum(0)
            .xyz()
            .cross(&event.outgoing_momentum(1).xyz());

        // Compute the cosine of the angle between the beam and this vector
        let cos_num = p_el.xyz().dot(&n_ppp);
        let cos_denom = p_el[E] * n_ppp.norm();

        // Check if the (beam, normal to photon plane) angle passes the cut
        if abs(cos_num) < self.beam_photon_plane_cut * cos_denom {
            return false;
        }

        // If all checks passed, we're good
        true
    }
}
