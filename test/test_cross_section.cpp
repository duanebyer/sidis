#include <catch2/catch.hpp>

#include <fstream>
#include <istream>
#include <sstream>
#include <utility>

#include <sidis/sidis.hpp>
#include <sidis/sf_model/ww.hpp>

#include "file_generator.hpp"
#include "rel_matcher.hpp"

using namespace sidis;

namespace {

struct Input {
	char particle_id;
	Real beam_energy;
	Real beam_pol;
	math::Vec3 target_pol;
	kin::PhaseSpace phase_space;
};

struct Output {
	Real born;
	Real amm;
	Real delta_vr;
	Real delta_vac_lep;
	Real delta_vac_had;
	Real err_born;
	Real err_amm;
	Real err_delta_vr;
	Real err_delta_vac_lep;
	Real err_delta_vac_had;
};

struct TestPair {
	Input input;
	Output output;
};

std::istream& operator>>(std::istream& in, TestPair& pair) {
	Input& input = pair.input;
	Output& output = pair.output;
	in >> input.particle_id;
	if (!(input.particle_id == 'e'
			|| input.particle_id == 'm'
			|| input.particle_id == 't')) {
		in.setstate(std::ios_base::failbit);
	}
	in >> input.beam_energy;
	in >> input.beam_pol;
	in >> input.target_pol.x;
	in >> input.target_pol.y;
	in >> input.target_pol.z;
	in >> input.phase_space.x;
	in >> input.phase_space.y;
	in >> input.phase_space.z;
	in >> input.phase_space.ph_t_sq;
	in >> input.phase_space.phi_h;
	in >> input.phase_space.phi;
	in >> output.born;
	in >> output.err_born;
	in >> output.amm;
	in >> output.err_amm;
	in >> output.delta_vr;
	in >> output.err_delta_vr;
	in >> output.delta_vac_lep;
	in >> output.err_delta_vac_lep;
	in >> output.delta_vac_had;
	in >> output.err_delta_vac_had;
	return in;
}

}

TEST_CASE(
		"Non-radiative cross-section values",
		"[xs]") {
	// Load the structure function data once for all tests.
	static sf::model::WW ww;

	// Load pre-computed data to use for cross-section verifications.
	TestPair test_pair = GENERATE(from_stream<TestPair>(std::move(
		std::ifstream("data/nrad_xs_vals.dat")), true));
	Input input = test_pair.input;
	Output output = test_pair.output;

	// Set up the input to the cross-section calculation.
	Real E_b = input.beam_energy;
	Real M = constant::MASS_P;
	Real m = 0.;
	if (input.particle_id == 'e') {
		m = constant::MASS_E;
	} else if (input.particle_id == 'm') {
		m = constant::MASS_MU;
	} else if (input.particle_id == 't') {
		m = constant::MASS_TAU;
	}
	Real mh = constant::MASS_PI;
	Real M_th = constant::MASS_P + constant::MASS_PI_0;
	Real pi = constant::PI;
	kin::Initial initial_state(M, m, E_b);
	kin::PhaseSpace phase_space = input.phase_space;
	kin::Kinematics kin(initial_state, phase_space, mh, M_th);
	// Get the final state particles.
	kin::Final final_state(initial_state, kin);
	// Get beam and target polarizations.
	Real beam_pol = input.beam_pol;
	math::Vec3 target_pol = input.target_pol;
	// Calculate structure functions.
	sf::Sf sf = ww.sf(kin.x, kin.z, kin.Q_sq, kin.ph_t);
	// Compute the cross-sections.
	Real born = xs::born(beam_pol, target_pol, kin, sf);
	Real amm = xs::amm(beam_pol, target_pol, kin, sf);
	Real delta_vr = xs::delta_vr(kin);
	Real delta_vac_lep = xs::delta_vac_lep(kin);
	Real delta_vac_had = xs::delta_vac_had(kin);

	// Print state information.
	std::stringstream ss;
	ss
		<< "pid   = " << input.particle_id   << std::endl
		<< "x     = " << phase_space.x       << std::endl
		<< "y     = " << phase_space.y       << std::endl
		<< "z     = " << phase_space.z       << std::endl
		<< "ph_t² = " << phase_space.ph_t_sq << std::endl
		<< "φ_h   = " << phase_space.phi_h   << std::endl
		<< "φ     = " << phase_space.phi     << std::endl
		<< "λ_e   = " << beam_pol            << std::endl
		<< "η_1   = " << target_pol.x        << std::endl
		<< "η_2   = " << target_pol.y        << std::endl
		<< "η_3   = " << target_pol.z;
	INFO(ss.str());

	// Do comparisons.
	CHECK_THAT(
		born,
		RelMatcher<Real>(output.born, output.err_born));
	CHECK_THAT(
		amm,
		RelMatcher<Real>(output.amm, output.err_amm));
	CHECK_THAT(
		delta_vr,
		RelMatcher<Real>(output.delta_vr, output.err_delta_vr));
	CHECK_THAT(
		delta_vac_lep,
		RelMatcher<Real>(output.delta_vac_lep, output.err_delta_vac_lep));
	CHECK_THAT(
		delta_vac_had,
		RelMatcher<Real>(output.delta_vac_had, output.err_delta_vac_had));
}

