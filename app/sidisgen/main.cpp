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
	kin::Initial initial_state {
		math::Vec4(M, 0., 0., 0.),
		math::Vec4(E_b, std::sqrt(E_b * E_b - m * m) * math::Vec3::Z),
	};
	// Choose the region of phase space to look at.
	kin::PhaseSpace phase_space {
		0.2,
		0.9,
		0.2,
		2.,
		5. * pi / 6.,
		-2. * pi / 6.
	};
	// Do the kinematic calculations.
	kin::Kinematics kin(initial_state, phase_space, mh, M_th);
	// Get the final state particles.
	kin::Final final_state(initial_state, kin);
	// Choose beam and target polarizations.
	Real beam_pol = 1.;
	math::Vec3 target_pol(0.0, 0.0, 1.0);
	// Transform the target polarization into a different basis.
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
	// Compute the cross-sections.
	Real born_cross_section = xs::born(beam_pol, eta, kin, sf);
	Real amm_cross_section = xs::amm(beam_pol, eta, kin, sf);
	std::cout << "Born cross-section is: " << born_cross_section << std::endl;
	std::cout << "AMM cross-section is:  " << amm_cross_section << std::endl;
	return 0;
}

