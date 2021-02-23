#include <catch2/catch.hpp>

#include <fstream>
#include <istream>
#include <memory>
#include <sstream>
#include <utility>

#include <sidis/sidis.hpp>
#include <sidis/sf_set/prokudin.hpp>
#include <sidis/sf_set/test.hpp>

#include "rel_matcher.hpp"
#include "stream_generator.hpp"

using namespace sidis;

namespace {

struct Input {
	int sf_set_idx;
	Real k0_cut;
	int beam_id;
	Real S;
	Real beam_pol;
	math::Vec3 target_pol;
	kin::PhaseSpace ph_space;
};

struct Output {
	Real born;
	Real err_born;
	Real amm;
	Real err_amm;
	Real nrad;
	Real err_nrad;
};

struct TestPair {
	Input input;
	Output output;
};

struct InputRad {
	int sf_set_idx;
	Real k0_cut;
	int beam_id;
	Real S;
	Real beam_pol;
	math::Vec3 target_pol;
	kin::PhaseSpaceRad ph_space;
};

struct OutputRad {
	Real rad_f;
	Real err_rad_f;
	Real rad;
	Real err_rad;
};

struct TestPairRad {
	InputRad input;
	OutputRad output;
};

std::istream& operator>>(std::istream& in, TestPair& pair) {
	Input& input = pair.input;
	Output& output = pair.output;
	in >> input.sf_set_idx;
	if (!(input.sf_set_idx >= -18 && input.sf_set_idx < 1)) {
		in.setstate(std::ios_base::failbit);
	}
	in >> input.k0_cut;
	in >> input.beam_id;
	if (!(input.beam_id >= 0 && input.beam_id < 3)) {
		in.setstate(std::ios_base::failbit);
	}
	in >> input.S;
	in >> input.beam_pol;
	input.target_pol.x = 0.;
	in >> input.target_pol.y;
	in >> input.target_pol.z;
	in >> input.ph_space.x;
	in >> input.ph_space.y;
	in >> input.ph_space.z;
	in >> input.ph_space.ph_t_sq;
	in >> input.ph_space.phi_h;
	in >> input.ph_space.phi;
	in >> output.born;
	in >> output.err_born;
	in >> output.amm;
	in >> output.err_amm;
	in >> output.nrad;
	in >> output.err_nrad;
	return in;
}

std::istream& operator>>(std::istream& in, TestPairRad& pair) {
	InputRad& input = pair.input;
	OutputRad& output = pair.output;
	in >> input.sf_set_idx;
	if (!(input.sf_set_idx >= -18 && input.sf_set_idx < 1)) {
		in.setstate(std::ios_base::failbit);
	}
	in >> input.k0_cut;
	in >> input.beam_id;
	if (!(input.beam_id >= 0 && input.beam_id < 3)) {
		in.setstate(std::ios_base::failbit);
	}
	in >> input.S;
	in >> input.beam_pol;
	input.target_pol.x = 0.;
	in >> input.target_pol.y;
	in >> input.target_pol.z;
	in >> input.ph_space.x;
	in >> input.ph_space.y;
	in >> input.ph_space.z;
	in >> input.ph_space.ph_t_sq;
	in >> input.ph_space.phi_h;
	in >> input.ph_space.phi;
	in >> input.ph_space.tau;
	in >> input.ph_space.phi_k;
	in >> input.ph_space.R;
	in >> output.rad_f;
	in >> output.err_rad_f;
	in >> output.rad;
	in >> output.err_rad;
	return in;
}

}

TEST_CASE(
		"Non-radiative cross-section values",
		"[xs]") {
	// Load pre-computed data to use for cross-section verifications.
	TestPair test_pair = GENERATE(
		from_stream<TestPair>(
			std::move(std::ifstream("data/xs_nrad_vals.dat")),
			true));
	Input input = test_pair.input;
	Output output = test_pair.output;

	std::unique_ptr<sf::SfSet> sf;
	if (input.sf_set_idx == 0) {
		sf.reset(new sf::set::ProkudinSfSet());
	} else {
		bool mask[18] = { false };
		mask[-input.sf_set_idx - 1] = true;
		sf.reset(new sf::set::TestSfSet(part::Nucleus::P, mask));
	}

	// Set up the input to the cross-section calculation.
	part::Lepton lep;
	if (input.beam_id == 0) {
		lep = part::Lepton::E;
	} else if (input.beam_id == 1) {
		lep = part::Lepton::MU;
	} else if (input.beam_id == 2) {
		lep = part::Lepton::TAU;
	}
	Real Mth = MASS_P + MASS_PI_0;
	kin::Particles ps(part::Nucleus::P, lep, part::Hadron::PI_P, Mth);
	kin::PhaseSpace ph_space = input.ph_space;
	kin::Kinematics kin(ps, input.S, ph_space);
	// Get beam and target polarizations.
	Real beam_pol = input.beam_pol;
	math::Vec3 eta = frame::hadron_from_target(kin) * input.target_pol;
	// Compute the cross-sections.
	Real born = xs::born(beam_pol, eta, kin, *sf);
	Real amm = xs::amm(beam_pol, eta, kin, *sf);
	Real nrad = xs::nrad_ir(beam_pol, eta, kin, *sf, input.k0_cut);

	// Print state information.
	std::stringstream ss;
	ss
		<< "sf_set_idx = " << input.sf_set_idx    << std::endl
		<< "pid        = " << part::name(lep) << std::endl
		<< "S          = " << input.S             << std::endl
		<< "x          = " << ph_space.x       << std::endl
		<< "y          = " << ph_space.y       << std::endl
		<< "z          = " << ph_space.z       << std::endl
		<< "ph_t²      = " << ph_space.ph_t_sq << std::endl
		<< "φ_h        = " << ph_space.phi_h   << std::endl
		<< "φ          = " << ph_space.phi     << std::endl
		<< "λ_e        = " << beam_pol            << std::endl
		<< "η_1        = " << eta.x               << std::endl
		<< "η_2        = " << eta.y               << std::endl
		<< "η_3        = " << eta.z;
	INFO(ss.str());

	// Do comparisons.
	CHECK_THAT(
		born,
		RelMatcher<Real>(output.born, 10.*output.err_born));
	CHECK_THAT(
		amm,
		RelMatcher<Real>(output.amm, 10.*output.err_amm));
	CHECK_THAT(
		nrad,
		RelMatcher<Real>(output.nrad, 10.*output.err_nrad));
}

TEST_CASE(
		"Radiative cross-section values",
		"[xs]") {
	// Load pre-computed data to use for cross-section verifications.
	TestPairRad test_pair = GENERATE(
		from_stream<TestPairRad>(
			std::move(std::ifstream("data/xs_rad_vals.dat")),
			true));
	InputRad input = test_pair.input;
	OutputRad output = test_pair.output;

	std::unique_ptr<sf::SfSet> sf;
	if (input.sf_set_idx == 0) {
		sf.reset(new sf::set::ProkudinSfSet());
	} else {
		bool mask[18] = { false };
		mask[-input.sf_set_idx - 1] = true;
		sf.reset(new sf::set::TestSfSet(part::Nucleus::P, mask));
	}

	// Set up the input to the cross-section calculation.
	part::Lepton lep;
	if (input.beam_id == 0) {
		lep = part::Lepton::E;
	} else if (input.beam_id == 1) {
		lep = part::Lepton::MU;
	} else if (input.beam_id == 2) {
		lep = part::Lepton::TAU;
	}
	Real Mth = MASS_P + MASS_PI_0;
	kin::Particles ps(part::Nucleus::P, lep, part::Hadron::PI_P, Mth);
	kin::PhaseSpaceRad ph_space = input.ph_space;
	kin::KinematicsRad kin(ps, input.S, ph_space);
	// Get beam and target polarizations.
	Real beam_pol = input.beam_pol;
	math::Vec3 eta = frame::hadron_from_target(kin.project()) * input.target_pol;
	// Compute the cross-sections.
	Real rad_f = xs::rad_f(beam_pol, eta, kin, *sf);
	Real rad = xs::rad(beam_pol, eta, kin, *sf);

	// Print state information.
	std::stringstream ss;
	ss
		<< "sf_set_idx = " << input.sf_set_idx    << std::endl
		<< "pid        = " << part::name(lep) << std::endl
		<< "S          = " << input.S             << std::endl
		<< "x          = " << ph_space.x       << std::endl
		<< "y          = " << ph_space.y       << std::endl
		<< "z          = " << ph_space.z       << std::endl
		<< "ph_t²      = " << ph_space.ph_t_sq << std::endl
		<< "φ_h        = " << ph_space.phi_h   << std::endl
		<< "φ          = " << ph_space.phi     << std::endl
		<< "λ_e        = " << beam_pol            << std::endl
		<< "η_1        = " << eta.x               << std::endl
		<< "η_2        = " << eta.y               << std::endl
		<< "η_3        = " << eta.z;
	INFO(ss.str());

	// Do comparisons.
	CHECK_THAT(
		rad_f,
		RelMatcher<Real>(output.rad_f, 10.*output.err_rad_f));
	CHECK_THAT(
		rad,
		RelMatcher<Real>(output.rad, 10.*output.err_rad));
}

