#include <iostream>
#include <sidis/sidis.hpp>
#include <sidis/sf_set/prokudin.hpp>

sidis::Real const PI = sidis::PI;
sidis::Real const M_TH = sidis::MASS_P + sidis::MASS_PI_0;

int main() {
	sidis::part::Particles particles(
		sidis::part::Nucleus::P,   // Target nucleus.
		sidis::part::Lepton::E,    // Beam lepton.
		sidis::part::Hadron::PI_P, // Leading hadron.
		M_TH                       // Threshold mass of undetected part.
	);
	sidis::Real S = 2. * 10.6 * particles.M; // Kinematic variable `S = 2 p k1`.
	sidis::kin::PhaseSpace phase_space {
		0.2,      // Bjorken x.
		0.9,      // Bjorken y.
		0.3,      // Bjorken z.
		2.,       // Transverse momentum of hadron, squared.
		0.5 * PI, // Azimuthal angle of hadron.
		0.,       // Azimuthal angle of transverse target polarization.
	};
	sidis::kin::Kinematics kin(particles, S, phase_space);
	sidis::Real beam_pol = 0.;
	sidis::math::Vec3 target_pol(0., 0., 0.);
	// Compute structure functions with WW-type approximation.
	sidis::sf::set::ProkudinSfSet sf;
	sidis::Real born_xs = sidis::xs::born(kin, sf, beam_pol, target_pol);
	std::cout << "Born unpolarized cross-section is " << born_xs << std::endl;
	return 0;
}

