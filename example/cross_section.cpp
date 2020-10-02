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
	Real M = MASS_P;
	Real m = MASS_E;
	Real mh = MASS_PI;
	Real M_th = MASS_P + MASS_PI;
	Initial initial_state(M, m, beam_energy);
	PhaseSpace phase_space { x, y, z, ph_t * ph_t, phi_h, phi };
	// Do kinematics calculations.
	Kinematics kin(initial_state, phase_space, mh, M_th);
	// Get the final state particles.
	Final final_state(initial_state, target_pol, kin);

	// Compute structure functions.
	sf::model::WW ww;
	sf::Sf sf = ww.sf(kin.x, kin.z, kin.Q_sq, kin.ph_t_sq);

	// Get the target polarization in the hadron frame.
	Vec3 eta = (frame::hadron_from_target(kin) * math::Vec4(0., target_pol)).r();
	// Compute cross-sections.
	Real born = xs::born(beam_pol, eta, kin, sf);
	Real born_rad_factor = xs::born_rad_factor(kin);
	Real amm = xs::amm(beam_pol, eta, kin, sf);
	Real tau_min = kin::KinematicsRad(kin, 0., 0., 0.).tau_min;
	Real tau_max = kin::KinematicsRad(kin, 0., 0., 0.).tau_max;
	// TODO: Fix the radiative cross-section calculation so it integrates more
	// precisely around the ridge in `(tau, phi_k)`.
	Real rad = integ::trapezoid([&](Real tau) {
			Real R_max = kin::KinematicsRad(kin, tau, 0., 0.).R_max;
			// Right now, the cross-section calculation breaks at `R=0`, so we
			// need to put in an arbitrary lower bound.
			// TODO: Remove the lower bound once cross-section calculation is
			// fixed.
			Real R_min = 0.001 * R_max;
			return integ::riemann(
				[&](Real phi_k) {
					return integ::trapezoid([&](Real R) {
							kin::KinematicsRad kin_rad(kin, tau, phi_k, R);
							sf::Sf sf_shift = ww.sf(
								kin_rad.shift_x,
								kin_rad.shift_z,
								kin_rad.shift_Q_sq,
								kin_rad.shift_ph_t_sq);
							Real rad = xs::rad(beam_pol, eta, kin_rad, sf, sf_shift);
							if (!std::isfinite(rad)) {
								// On occasion the shifted kinematics take us
								// out of the region on which grid data exists
								// for the structure functions. For now, just
								// set to zero when that happens.
								// TODO: Deal with the boundary properly.
								return 0.;
							} else {
								return rad;
							}
						},
						R_min, R_max, 50);
				},
				0., 2. * PI, 4);
		},
		tau_min, tau_max, 50);

	std::cout << "σ_B   = " << born << std::endl;
	std::cout << "σ_AMM = " << amm << std::endl;
	std::cout << "σ_rad = " << rad << std::endl;
	std::cout << "σ_tot = " << born_rad_factor * born + amm + rad << std::endl;

	return 0;
}

