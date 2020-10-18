#include <cmath>
#include <exception>
#include <iostream>
#include <stdexcept>
#include <string>

#include <sidis/sidis.hpp>
#include <sidis/sf_model/ww.hpp>
#include <sidis/extra/integrate.hpp>
#include <sidis/extra/transform.hpp>
#include <sidis/extra/vector.hpp>

using namespace sidis;
using namespace sidis::constant;
using namespace sidis::kin;
using namespace sidis::math;

int main(int argc, char** argv) {
	// Read input parameters from command line.
	Real beam_energy;
	Real beam_pol;
	Vec3 target_pol;
	Real x, y, z, ph_t, phi_h, phi;
	try {
		if (argc != 10) {
			throw std::invalid_argument("Expected 9 command line arguments");
		}
		beam_energy = std::stold(argv[1]);
		std::string beam_pol_str = std::string(argv[2]);
		std::string target_pol_str = std::string(argv[3]);
		x = std::stold(argv[4]);
		y = std::stold(argv[5]);
		z = std::stold(argv[6]);
		ph_t = std::stold(argv[7]);
		phi_h = std::stold(argv[8]);
		phi = std::stold(argv[9]);
		if (beam_pol_str == "U") {
			beam_pol = 0.;
		} else if (beam_pol_str == "L") {
			beam_pol = 1.;
		} else {
			throw std::out_of_range(
				"Beam must be unpolarized (U) or longitudinally polarized (L)");
		}
		if (target_pol_str == "U") {
			target_pol = Vec3::ZERO;
		} else if (target_pol_str == "L") {
			target_pol = Vec3::Z;
		} else if (target_pol_str == "T") {
			target_pol = Vec3::Y;
		} else {
			throw std::out_of_range(
				"Target must be unpolarized (U), longitudinally polarized (L), "
				"or tangentially polarized (T)");
		}
	} catch (std::exception const& e) {
		std::cerr << "Error: " << e.what() << std::endl;
		std::cout << "Usage: "
			<< "cross-section "
			<< "<beam energy> "
			<< "<beam pol. U,L> "
			<< "<target pol. U,L,T> "
			<< "<x> <y> <z> <ph_t> <φ_h> <φ>"
			<< std::endl;
		return 1;
	}

	// Set up the initial state particles.
	Real M_th = MASS_P + MASS_PI;
	Initial initial_state(Nucleus::P, Lepton::E, beam_energy);
	PhaseSpace phase_space { x, y, z, ph_t * ph_t, phi_h, phi };
	// Do kinematics calculations.
	Kinematics kin(initial_state, phase_space, Hadron::PI_P, M_th);
	// Get the final state particles.
	Final final_state(initial_state, target_pol, kin);

	// Construct a model for computing structure functions.
	sf::model::WW ww;

	// Get the target polarization in the hadron frame.
	Vec3 eta = frame::hadron_from_target(kin) * target_pol;
	// Compute cross-sections.
	Real born = xs::born(beam_pol, eta, kin, ww);
	Real born_rad_factor = xs::born_rad_factor(kin);
	Real amm = xs::amm(beam_pol, eta, kin, ww);
	Real rad = xs::rad_integ(beam_pol, eta, kin, ww);
	std::cout << "σ_B   = " << born << std::endl;
	std::cout << "σ_AMM = " << amm << std::endl;
	std::cout << "σ_rad = " << rad << std::endl;
	std::cout << "σ_tot = " << born_rad_factor * born + amm + rad << std::endl;

	return 0;
}

