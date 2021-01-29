#include <iostream>
#include <sidis/sidis.hpp>
#include <sidis/sf_set/ww.hpp>
#include <sidis/extra/vector.hpp>

sidis::Real const PI = sidis::constant::PI;
sidis::Real const M_TH = sidis::constant::MASS_P + sidis::constant::MASS_PI_0;

int main() {
	sidis::kin::Particles particles(
		sidis::constant::Nucleus::P,   // Target nucleus.
		sidis::constant::Lepton::E,    // Beam lepton.
		sidis::constant::Hadron::PI_P, // Leading hadron.
		M_TH                           // Threshold mass of undetected part.
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
	sidis::sf::model::WW ww;
	sidis::Real born_xs = sidis::xs::born(beam_pol, target_pol, kin, ww);
	std::cout << "Born unpolarized cross-section is " << born_xs << std::endl;
	return 0;
}

