//! Mechanism to apply a cut to generated events

use ::{
    event::{Event, OUTGOING_COUNT},
    linalg::{
        dimension::*,
        momentum::{E, X, xyz},
    },
    numeric::{
        functions::*,
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

        // Check if the (photon1, photon{2, 3}) angles pass the cut
        // FIXME: Turn this back into a loop once const generics allow for it
        assert_eq!(OUTGOING_COUNT, 3, "This part assumes 3 outgoing particles");
        let p_ph1 = event.outgoing_momentum(0);
        let p_ph23 = p_out.fixed_rows::<U2>(1);
        let p_ph23_xyz = p_ph23.fixed_columns::<U3>(X);
        let p_ph23_e = p_ph23.fixed_columns::<U1>(E);
        let cos_num_1x23 = p_ph23_xyz * xyz(&p_ph1);
        let cos_denom_1x23 = p_ph23_e * p_ph1[E];
        for (&num, denom) in cos_num_1x23.iter().zip(cos_denom_1x23.iter()) {
            if num > self.b_cut * denom { return false; }
        }

        // Check if the (photon2, photon3) angle passes the cut
        // FIXME: Merge with the above loop once we can have it
        let p_ph2 = p_ph23.fixed_rows::<U1>(0);
        let p_ph3 = p_ph23.fixed_rows::<U1>(1);
        let p_ph2_xyz = p_ph2.fixed_columns::<U3>(X);
        let p_ph3_xyz = p_ph3.fixed_columns::<U3>(X);
        let cos_num_2x3 = p_ph2_xyz.dot(&p_ph3_xyz);
        let cos_denom_2x3 = p_ph2[E] * p_ph3[E];
        if cos_num_2x3 > self.b_cut * cos_denom_2x3 { return false; }

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
