#include <cmath>
#include <iostream>

#include <algorithm>
#include <sidis/sidis.hpp>
#include <sidis/sf_model/ww.hpp>

using namespace sidis;

int main(int argc, char** argv) {
	// Set up the initial state.
	Real E_b = 100;
	Real M = constant::MASS_P;
	Real m = constant::MASS_E;
	Real mh = constant::MASS_PI;
	Real M_th = constant::MASS_P + constant::MASS_PI;
	Real pi = constant::PI;
	kin::Initial initial_state(M, m, E_b);
	// Choose the region of phase space to look at.
	kin::PhaseSpace phase_space {
		0.2,
		0.9,
		0.2,
		2.,
		5. * pi / 6.,
		-2. * pi / 6.,
	};
	// Do the kinematic calculations.
	kin::Kinematics kin(initial_state, phase_space, mh, M_th);
	kin::KinematicsRad kin_rad(kin, 0.2, 2. * pi / 6., 1.);
	// Get the final state particles.
	kin::Final final_state(initial_state, kin);
	kin::FinalRad final_state_rad(initial_state, kin_rad);
	// Choose beam and target polarizations.
	Real beam_pol = 1.;
	math::Vec3 target_pol(0.0, 0.0, 1.0);
	// Transform the target polarization into a different basis (same for
	// radiative and non-radiative parts).
	math::Vec3 e_z = final_state.q.r.unit();
	math::Vec3 e_y = math::cross(e_z, final_state.ph.r).unit();
	math::Vec3 e_x = math::cross(e_y, e_z).unit();
	math::Vec3 eta(
		math::dot(target_pol, e_x),
		math::dot(target_pol, e_y),
		math::dot(target_pol, e_z));
	// Choose structure functions. For now, use the WW approximation.
	sf::model::WW ww;
	sf::Sf sf = ww.sf(kin.x, kin.z, kin.Q_sq, kin.ph_t);
	sf::Sf shift_sf = ww.sf(kin_rad.shift_x, kin_rad.shift_z, kin_rad.shift_Q_sq, kin_rad.shift_ph_t);
	// Compute the cross-sections.
	Real born_xs = xs::born(beam_pol, eta, kin, sf);
	Real amm_xs = xs::amm(beam_pol, eta, kin, sf);
	Real rad_xs = xs::rad(beam_pol, eta, kin_rad, sf, shift_sf);
	std::cout << "Born cross-section is:      " << born_xs << std::endl;
	std::cout << "AMM cross-section is:       " << amm_xs << std::endl;
	std::cout << "Radiative cross-section is: " << rad_xs << std::endl;
	return 0;
}

