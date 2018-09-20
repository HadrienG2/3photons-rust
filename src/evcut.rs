//! Mechanism to apply a cut to generated events

use ::{
    event::{Event, OUTGOING_COUNT},
    linalg::{E, U1, U3, X, xyz},
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
        let p_out_xyz = p_out.fixed_columns::<U3>(X);
        let p_out_e = p_out.fixed_columns::<U1>(E);

        // Check if the (beam, photon) angles pass the cut
        {
            let cos_num = p_out_xyz * xyz(&p_el);
            let cos_denom = p_out_e * p_el[E];
            for (&num, denom) in cos_num.iter().zip(cos_denom.iter()) {
                if abs(num) > self.a_cut * denom { return false; }
            }
        }

        // Check if the (photon, photon) angles pass the cut
        for par1 in 0..OUTGOING_COUNT-1 {
            // Pick one photon
            let p_ph = event.outgoing_momentum(par1);

            // Pick the other photons coming after it (we want to study pairs)
            let other_p_ph = p_out.rows(par1+1, OUTGOING_COUNT-par1-1);
            let other_p_ph_xyz = other_p_ph.fixed_columns::<U3>(X);
            let other_p_ph_e = other_p_ph.fixed_columns::<U1>(E);

            // Do the same thing as for (beam, photon) angles
            let cos_num = other_p_ph_xyz * xyz(&p_ph);
            let cos_denom = other_p_ph_e * p_ph[E];
            for (&num, denom) in cos_num.iter().zip(cos_denom.iter()) {
                if abs(num) > self.b_cut * denom { return false; }
            }
        }

        // Compute a vector which is normal to the outgoing photon plane
        // NOTE: This notion is only valid because we have three output photons
        debug_assert_eq!(OUTGOING_COUNT, 3);
        let n_ppp = xyz(&event.outgoing_momentum(0))
                        .cross(&xyz(&event.outgoing_momentum(1)));

        // Compute the cosine of the angle between the beam and this vector
        let cos_num = xyz(&p_el).dot(&n_ppp);
        let cos_denom = p_el[E] * n_ppp.norm();

        // Check if the (beam, normal to photon plane) angle passes the cut
        abs(cos_num) >= self.sin_cut * cos_denom
    }
}
