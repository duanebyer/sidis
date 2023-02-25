#include <iomanip>
#include <ios>
#include <iostream>
#include <memory>
#include <stdexcept>
#include <string>

#include <sidis/sidis.hpp>
#include <sidis/sf_set/prokudin.hpp>
#include <sidis/sf_set/mask.hpp>
#include <sidis/sf_set/test.hpp>

using namespace sidis;
using namespace sidis::kin;
using namespace sidis::math;
using namespace sidis::part;

// This program calculates the cross-section given the parameters passed to it
// on the command line. It outputs the Born cross-section, and the various
// radiative correction contributions.
int main(int argc, char** argv) {
	Real Mth = MASS_P + MASS_PI_0;
	part::Lepton beam = part::Lepton::E;
	part::Nucleus target = part::Nucleus::P;
	part::Hadron hadron = part::Hadron::PI_P;

	// Read input parameters from command line.
	Real beam_energy;
	Real beam_pol;
	Vec3 target_pol;
	std::unique_ptr<sidis::sf::SfSet> sf;
	std::unique_ptr<sidis::sf::TmdSet> tmd;
	Real x, y, z, ph_t_sq, phi_h, phi;
	Real k0_cut = 0.01;
	Real tau = 0., phi_k = 0., R = 0.;
	bool radiative;
	math::IntegParams params { math::IntegMethod::CUBATURE, 100000, 10000 };
	try {
		if (argc != 12 && argc != 14) {
			throw std::invalid_argument(
				"Unexpected number of command line arguments");
		}
		int set_idx = std::stoi(argv[1]);
		beam_energy = std::stold(argv[2]);
		std::string beam_pol_str = std::string(argv[3]);
		std::string target_pol_str = std::string(argv[4]);
		x = std::stold(argv[5]);
		y = std::stold(argv[6]);
		z = std::stold(argv[7]);
		ph_t_sq = std::stold(argv[8]);
		phi_h = std::stold(argv[9]);
		phi = std::stold(argv[10]);
		if (argc == 12) {
			radiative = false;
			k0_cut = std::stold(argv[11]);
		} else {
			radiative = true;
			tau = std::stold(argv[11]);
			phi_k = std::stold(argv[12]);
			R = std::stold(argv[13]);
		}

		if (set_idx == 0) {
			sf.reset(new sf::set::ProkudinSfSet());
		} else if (set_idx <= -1 && set_idx >= -static_cast<int>(sf::set::NUM_SF)) {
			bool mask[sf::set::NUM_SF] = { false };
			mask[-set_idx - 1] = true;
			sf.reset(new sf::set::MaskSfSet(
				mask, new sf::set::TestSfSet(target)));
		} else {
			throw std::out_of_range(
				"SF set index must be Prokudin (0), STF (1), or Test ("
				+ std::to_string(-sf::set::NUM_SF) + " to -1)");
		}

		if (beam_pol_str == "U") {
			beam_pol = 0.;
		} else if (beam_pol_str == "L") {
			beam_pol = 1.;
		} else {
			throw std::out_of_range(
				"Beam must be unpolarized (U) or longitudinally polarized (L)");
		}
		if (target_pol_str == "U") {
			target_pol = VEC3_ZERO;
		} else if (target_pol_str == "L") {
			target_pol = VEC3_Z;
		} else if (target_pol_str == "T") {
			target_pol = VEC3_Y;
		} else {
			throw std::out_of_range(
				"Target must be unpolarized (U), longitudinally polarized (L), "
				"or tangentially polarized (T)");
		}
	} catch (std::exception const& e) {
		std::cerr << "Error: " << e.what() << std::endl;
		std::cout << "Usage" << std::endl
			<< "non-radiative: "
			<< "cross_section "
			<< "<SF set idx> "
			<< "<E_b> "
			<< "<U,L> "
			<< "<U,L,T> "
			<< "<x> <y> <z> <ph_t²> <φ_h> <φ> "
			<< "<k0 cut>" << std::endl
			<< "radiative:     "
			<< "cross_section "
			<< "<SF set idx> "
			<< "<E_b> "
			<< "<U,L> "
			<< "<U,L,T> "
			<< "<x> <y> <z> <ph_t²> <φ_h> <φ> "
			<< "<τ> <φ_k> <R>" << std::endl;
		return 1;
	}

	// Set up the initial state particles.
	Particles ps(target, beam, hadron, Mth);
	Real S = 2. * ps.M * beam_energy;
	PhaseSpace ph_space { x, y, z, ph_t_sq, phi_h, phi };
	// Check that we are in valid kinematic region.
	if (!(S >= cut::S_min(ps))) {
		throw std::out_of_range("Beam energy is below threshold");
	} else if (!cut::x_bound(ps, S).contains(x)) {
		throw std::out_of_range("x is out of valid kinematic range");
	} else if (!cut::y_bound(ps, S, x).contains(y)) {
		throw std::out_of_range("y is out of valid kinematic range");
	} else if (!cut::z_bound(ps, S, x, y).contains(z)) {
		throw std::out_of_range("z is out of valid kinematic range");
	} else if (!cut::ph_t_sq_bound(ps, S, x, y, z).contains(ph_t_sq)) {
		throw std::out_of_range("ph_t_sq is out of valid kinematic range");
	}
	// Do kinematics calculations.
	Kinematics kin(ps, S, ph_space);
	// Get the target polarization in the hadron frame.
	Vec3 eta = frame::hadron_from_target(kin) * target_pol;

	// Compute cross-sections.
	if (!radiative) {
		std::cout << std::scientific << std::setprecision(16);
		Real born = xs::born(kin, *sf, beam_pol, eta);
		std::cout << "σ_B       = " << born << std::endl;
		Real amm = xs::amm(kin, *sf, beam_pol, eta);
		std::cout << "σ_AMM     = " << amm << std::endl;
		Real nrad_ir = xs::nrad_ir(kin, *sf, beam_pol, eta, k0_cut);
		std::cout << "σ_nrad_ir = " << nrad_ir << std::endl;
		EstErr rad_f = xs::rad_f_integ(kin, *sf, beam_pol, eta, k0_cut, params);
		std::cout << "σ_rad_f   = " << rad_f.val << " ± " << rad_f.err << std::endl;
		EstErr nrad = xs::nrad_integ(kin, *sf, beam_pol, eta, k0_cut, params);
		std::cout << "σ_nrad    = " << nrad.val << " ± " << nrad.err << std::endl;
		EstErr rad = xs::rad_integ(kin, *sf, beam_pol, eta, k0_cut, params);
		std::cout << "σ_rad     = " << rad.val << " ± " << rad.err << std::endl;
		std::cout << "σ_tot     = " << nrad.val + rad.val << " ± " << nrad.err + rad.err << std::endl;
	} else {
		// Do radiative kinematics checks.
		if (!cut::tau_bound(kin).contains(tau)) {
			throw std::out_of_range("tau is out of valid kinematic range");
		} else if (!cut::R_bound(kin, tau, phi_k).contains(R)) {
			throw std::out_of_range("R is out of valid kinematic range");
		}
		std::cout << std::scientific << std::setprecision(16);
		KinematicsRad kin_rad(kin, tau, phi_k, R);
		Real rad_f = xs::rad_f(kin_rad, *sf, beam_pol, eta);
		std::cout << "σ_rad_f = " << rad_f << std::endl;
		Real rad = xs::rad(kin_rad, *sf, beam_pol, eta);
		std::cout << "σ_rad   = " << rad << std::endl;
	}

	return 0;
}

