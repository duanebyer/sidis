#include <iomanip>
#include <ios>
#include <iostream>
#include <random>
#include <stdexcept>
#include <string>

#include <sidis/sidis.hpp>

using namespace sidis;
using namespace sidis::kin;
using namespace sidis::math;
using namespace sidis::part;

// This program returns a random valid point in phase space.
int main(int argc, char** argv) {
	Real Mth = MASS_P + MASS_PI_0;
	part::Lepton beam = part::Lepton::E;
	part::Nucleus target = part::Nucleus::P;
	part::Hadron hadron = part::Hadron::PI_P;

	// Read input parameters from command line.
	Real beam_energy;
	bool radiative;
	try {
		if (argc != 3) {
			throw std::invalid_argument(
				"Unexpected number of command line arguments");
		}
		beam_energy = std::stold(argv[1]);
		std::string radiative_str = argv[2];
		if (radiative_str == "true"
				|| radiative_str == "on"
				|| radiative_str == "rad") {
			radiative = true;
		} else if (radiative_str == "false"
				|| radiative_str == "off"
				|| radiative_str == "nrad") {
			radiative = false;
		} else {
			throw std::out_of_range(
				"Must select radiative (rad) or "
				"non-radiative (nrad) phase space point");
		}
	} catch (std::exception const& e) {
		std::cerr << "Error: " << e.what() << std::endl;
		std::cout << "Usage: "
			<< "random_phase_space <E_b> <rad,nrad>"
			<< std::endl;
		return 1;
	}

	// Repeatedly choose a random point within phase space until we get one that
	// is kinematically valid.
	Particles ps(target, beam, hadron, Mth);
	Real S = 2. * ps.M * beam_energy;
	std::random_device rd;
	std::mt19937 rng(rd());
	std::uniform_real_distribution<Real> dist(0., 1.);

	Real point[9] = {
		dist(rng), dist(rng), dist(rng),
		dist(rng), dist(rng), dist(rng),
		dist(rng), dist(rng), dist(rng),
	};
	Kinematics kin;
	KinematicsRad kin_rad;
	while (!cut::take(ps, S, point, &kin, nullptr)) { }
	while (!cut::take(kin, point + 6, &kin_rad, nullptr)) { }

	std::cout << std::scientific << std::setprecision(16);

	std::cout << "x     = " << kin.x << std::endl;
	std::cout << "y     = " << kin.y << std::endl;
	std::cout << "z     = " << kin.z << std::endl;
	std::cout << "ph_t² = " << kin.ph_t_sq << std::endl;
	std::cout << "φ_h   = " << kin.phi_h << std::endl;
	std::cout << "φ     = " << kin.phi << std::endl;
	if (radiative) {
		std::cout << "τ     = " << kin_rad.tau << std::endl;
		std::cout << "φ_k   = " << kin_rad.phi_k << std::endl;
		std::cout << "R     = " << kin_rad.R << std::endl;
	}

	return 0;
}

