//! Mechanism to apply a cut to generated events

use ::{
    event::{Event, OUTGOING_COUNT},
    linalg::{self, E},
    numeric::{
        functions::abs,
        Real
    },
};


/// Cuts on generated events
pub struct EventCut {
    /// Cut on maximum cosine of (beam, photons) angle
    pub a_cut: Real,

    /// Cut on maximum cosine of (photon, photon) angle
    pub b_cut: Real,

    /// Cut on minimum photon energy
    pub e_min: Real,

    /// Cut on minimum cosine of (beam, normal to the photon plane) angle
    pub sin_cut: Real,
}
//
impl EventCut {
    /// Setup the cuts on generated events
    pub fn new(a_cut: Real, b_cut: Real, e_min: Real, sin_cut: Real) -> Self {
        EventCut {
            a_cut: a_cut,
            b_cut: b_cut,
            e_min: e_min,
            sin_cut: sin_cut,
        }
    }

    /// Decide whether a generated event passes the cut or should be rejected
    pub fn keep(&self, event: &Event) -> bool
    {
        // Check if the outgoing photons pass the energy cut
        if event.min_photon_energy() < self.e_min { return false; }

        // Get the incoming electron 4-momentum and outgoing photon 4-momenta
        let p_el = event.electron_momentum();
        let p_out = event.outgoing_momenta();

        // Compute the cosines of the angles between photons and the e- beam
        let mut cos_p_el = [0.; OUTGOING_COUNT];
        let p_el_xyz = linalg::xyz(&p_el);
        for (cos, p_ph) in cos_p_el.iter_mut().zip(p_out.iter()) {
            let p_ph_xyz = linalg::xyz(&p_ph);
            *cos = p_el_xyz.dot(&p_ph_xyz) / (p_el[E] * p_ph[E]);
        }

        // Check if the (beam, photon) angles pass the cut
        if cos_p_el.iter().any(|&cos| abs(cos) > self.a_cut) { return false; }

        // Compute the cosines of the angles between photon pairs
        const OUTGOING_PAIRS: usize = OUTGOING_COUNT * (OUTGOING_COUNT-1) / 2;
        let mut cos_p_p = [0.; OUTGOING_PAIRS];
        {
            let mut cos_iter = cos_p_p.iter_mut();
            for j in 1..OUTGOING_COUNT {
                for i in 0..j {
                    let cos = cos_iter.next().unwrap();
                    let p_ph_1 = &p_out[i];
                    let p_ph_2 = &p_out[j];
                    *cos = linalg::xyz(&p_ph_1).dot(&linalg::xyz(&p_ph_2)) /
                               (p_ph_1[E] * p_ph_2[E]);
                }
            }
        }

        // Check if the (photon, photon) angles pass the cut
        if cos_p_p.iter().any(|&cos| cos > self.b_cut) { return false; }

        // Compute a vector which is normal to the outgoing photon plane
        // NOTE: This notion is only valid because we have three output photons
        debug_assert_eq!(OUTGOING_COUNT, 3);
        let n_ppp = linalg::xyz(&p_out[0]).cross(&linalg::xyz(&p_out[1]));

        // Compute the cosine of the angle between the beam and this vector
        let cos_n = linalg::xyz(p_el).dot(&n_ppp) / (p_el[E] * n_ppp.norm());

        // Check if the (beam, normal to photon plane) angle passes the cut
        abs(cos_n) >= self.sin_cut
    }
}
