#include <cmath>
#include <iostream>

#include "sidis/sidis.hpp"

using namespace sidis;

int main(int argc, char** argv) {
	// Set up the initial state.
	Real E_b = 100;
	Real M = constant::MASS_P;
	Real m = constant::MASS_E;
	Real mh = constant::MASS_PI;
	Real M_th = constant::MASS_P + constant::MASS_PI;
	Real pi = constant::PI;
	kin::Initial initial {
		math::Vec4(M, 0., 0., 0.),
		math::Vec4(E_b, std::sqrt(E_b * E_b - m * m) * math::Vec3::Z),
	};

	// Choose the region of phase space to look at.
	kin::PhaseSpace phase_space {
		0.2,
		0.9,
		0.2,
		pi / 6.,
		5. * pi / 6.,
		-2. * pi / 6.
	};
	// Do the kinematic calculations.
	kin::Kinematics kin(initial, phase_space, mh, M_th);
	// Choose structure functions. For now, just use constants.
	sf::Sf sf(sf::SfUU { 0., 1., 0., 0. }, {}, {}, {}, {}, {});
	// Compute the cross-section.
	Real born_cross_section = xs::born(0., math::Vec3::ZERO, kin, sf);

	std::cout << "Born cross-section is: " << born_cross_section << std::endl;
	return 0;
}

