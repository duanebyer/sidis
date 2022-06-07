#include <catch2/catch.hpp>

#include <fstream>
#include <istream>
#include <limits>
#include <memory>
#include <sstream>
#include <utility>

#include <sidis/sidis.hpp>
#include <sidis/sf_set/mask.hpp>
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
	if (!(input.sf_set_idx >= -sf::set::NUM_SF && input.sf_set_idx < 1)) {
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
	if (!(input.sf_set_idx >= -sf::set::NUM_SF && input.sf_set_idx < 1)) {
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
		"Gaussian-WW structure functions",
		"[xs]") {
	// Verify that these two methods for computing structure functions are
	// equivalent.
	static sf::set::ProkudinTmdSet tmd_set;
	static sf::GaussianWwTmdSfSet sf_set_1(tmd_set);
	static sf::set::ProkudinSfSet sf_set_2;

	Real x = GENERATE(0.12, 0.31, 0.43, 0.23, 0.65, 0.74);
	Real Q_sq = GENERATE(2.4, 5.3, 7.9, 12.5);
	Real z = GENERATE(0.32, 0.54, 0.61, 0.68);
	Real ph_t = GENERATE(0.012, 0.102, 0.046, 0.242);
	Real ph_t_sq = ph_t * ph_t;

	part::Hadron h = part::Hadron::PI_P;

	sf::Sf sf_1 = sf_set_1.sf(h, x, z, Q_sq, ph_t_sq);
	sf::Sf sf_2 = sf_set_2.sf(h, x, z, Q_sq, ph_t_sq);

	Real prec = 1e2*std::numeric_limits<Real>::epsilon();
	// Two of the structure functions are commented out. For technical reasons,
	// those two are only approximately equal instead of exactly equal.
	CHECK_THAT(sf_1.F_UUT, RelMatcher<Real>(sf_2.F_UUT, prec));
	CHECK_THAT(sf_1.F_UU_cos_phih, RelMatcher<Real>(sf_2.F_UU_cos_phih, prec));
	CHECK_THAT(sf_1.F_UU_cos_2phih, RelMatcher<Real>(sf_2.F_UU_cos_2phih, prec));
	CHECK_THAT(sf_1.F_UL_sin_phih, RelMatcher<Real>(sf_2.F_UL_sin_phih, prec));
	CHECK_THAT(sf_1.F_UL_sin_2phih, RelMatcher<Real>(sf_2.F_UL_sin_2phih, prec));
	CHECK_THAT(sf_1.F_UTT_sin_phih_m_phis, RelMatcher<Real>(sf_2.F_UTT_sin_phih_m_phis, prec));
	//CHECK_THAT(sf_1.F_UT_sin_2phih_m_phis, RelMatcher<Real>(sf_2.F_UT_sin_2phih_m_phis, prec));
	CHECK_THAT(sf_1.F_UT_sin_3phih_m_phis, RelMatcher<Real>(sf_2.F_UT_sin_3phih_m_phis, prec));
	//CHECK_THAT(sf_1.F_UT_sin_phis, RelMatcher<Real>(sf_2.F_UT_sin_phis, prec));
	CHECK_THAT(sf_1.F_UT_sin_phih_p_phis, RelMatcher<Real>(sf_2.F_UT_sin_phih_p_phis, prec));
	CHECK_THAT(sf_1.F_LL, RelMatcher<Real>(sf_2.F_LL, prec));
	CHECK_THAT(sf_1.F_LL_cos_phih, RelMatcher<Real>(sf_2.F_LL_cos_phih, prec));
	CHECK_THAT(sf_1.F_LT_cos_phih_m_phis, RelMatcher<Real>(sf_2.F_LT_cos_phih_m_phis, prec));
	CHECK_THAT(sf_1.F_LT_cos_2phih_m_phis, RelMatcher<Real>(sf_2.F_LT_cos_2phih_m_phis, prec));
	CHECK_THAT(sf_1.F_LT_cos_phis, RelMatcher<Real>(sf_2.F_LT_cos_phis, prec));
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
		bool mask[sf::set::NUM_SF] = { false };
		mask[-input.sf_set_idx - 1] = true;
		sf.reset(new sf::set::MaskSfSet(
			mask, new sf::set::TestSfSet(part::Nucleus::P)));
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
	part::Particles ps(part::Nucleus::P, lep, part::Hadron::PI_P, Mth);
	kin::PhaseSpace ph_space = input.ph_space;
	kin::Kinematics kin(ps, input.S, ph_space);
	// Get beam and target polarizations.
	Real lambda = input.beam_pol;
	math::Vec3 eta = frame::hadron_from_target(kin) * input.target_pol;
	// Compute the cross-sections.
	// Note that we use a specific value for alpha here, which isn't the most
	// realistic (at Q^2=0, not evolved). It doesn't matter for these tests
	// though, as we simply want to protect against any changes to the
	// underlying cross-section math.
	Real alpha = 7.2973525664e-3L;
	xs::Born b_born(kin, alpha);
	xs::Amm b_amm(kin, alpha);
	xs::Nrad b_nrad(kin, input.k0_cut, alpha);
	had::HadXX had(kin, *sf);
	lep::LepBornXX lep_born(kin);
	lep::LepAmmXX lep_amm(kin);
	lep::LepNradXX lep_nrad(lep_born, lep_amm);

	Real born = xs::born_uu_base(b_born, lep_born, had)
		+ dot(eta, xs::born_up_base(b_born, lep_born, had))
		+ lambda * xs::born_lu_base(b_born, lep_born, had)
		+ lambda * dot(eta, xs::born_lp_base(b_born, lep_born, had));
	Real amm = xs::amm_uu_base(b_amm, lep_amm, had)
		+ dot(eta, xs::amm_up_base(b_amm, lep_amm, had))
		+ lambda * xs::amm_lu_base(b_amm, lep_amm, had)
		+ lambda * dot(eta, xs::amm_lp_base(b_amm, lep_amm, had));
	Real nrad = xs::nrad_ir_uu_base(b_nrad, lep_nrad, had)
		+ dot(eta, xs::nrad_ir_up_base(b_nrad, lep_nrad, had))
		+ lambda * xs::nrad_ir_lu_base(b_nrad, lep_nrad, had)
		+ lambda * dot(eta, xs::nrad_ir_lp_base(b_nrad, lep_nrad, had));

	// Print state information.
	std::stringstream ss;
	ss
		<< "sf_set_idx = " << input.sf_set_idx << std::endl
		<< "pid        = " << part::name(lep)  << std::endl
		<< "S          = " << input.S          << std::endl
		<< "x          = " << ph_space.x       << std::endl
		<< "y          = " << ph_space.y       << std::endl
		<< "z          = " << ph_space.z       << std::endl
		<< "ph_t²      = " << ph_space.ph_t_sq << std::endl
		<< "φ_h        = " << ph_space.phi_h   << std::endl
		<< "φ          = " << ph_space.phi     << std::endl
		<< "λ_e        = " << lambda           << std::endl
		<< "η_1        = " << eta.x            << std::endl
		<< "η_2        = " << eta.y            << std::endl
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
		bool mask[sf::set::NUM_SF] = { false };
		mask[-input.sf_set_idx - 1] = true;
		sf.reset(new sf::set::MaskSfSet(
			mask, new sf::set::TestSfSet(part::Nucleus::P)));
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
	part::Particles ps(part::Nucleus::P, lep, part::Hadron::PI_P, Mth);
	kin::PhaseSpaceRad ph_space = input.ph_space;
	kin::KinematicsRad kin(ps, input.S, ph_space);
	// Get beam and target polarizations.
	Real lambda = input.beam_pol;
	math::Vec3 eta = frame::hadron_from_target(kin.project()) * input.target_pol;
	// Compute the cross-sections.
	Real alpha = 7.2973525664e-3L;
	xs::Rad b_rad(kin, alpha);
	had::HadRadXX had_rad(kin, *sf);
	had::HadRadFXX had_rad_f(kin, *sf);
	lep::LepRadXX lep_rad(kin);

	Real rad = xs::rad_uu_base(b_rad, lep_rad, had_rad)
		+ dot(eta, xs::rad_up_base(b_rad, lep_rad, had_rad))
		+ lambda * xs::rad_lu_base(b_rad, lep_rad, had_rad)
		+ lambda * dot(eta, xs::rad_lp_base(b_rad, lep_rad, had_rad));
	Real rad_f = xs::rad_f_uu_base(b_rad, lep_rad, had_rad_f)
		+ dot(eta, xs::rad_f_up_base(b_rad, lep_rad, had_rad_f))
		+ lambda * xs::rad_f_lu_base(b_rad, lep_rad, had_rad_f)
		+ lambda * dot(eta, xs::rad_f_lp_base(b_rad, lep_rad, had_rad_f));

	// Print state information.
	std::stringstream ss;
	ss
		<< "sf_set_idx = " << input.sf_set_idx << std::endl
		<< "pid        = " << part::name(lep)  << std::endl
		<< "S          = " << input.S          << std::endl
		<< "x          = " << ph_space.x       << std::endl
		<< "y          = " << ph_space.y       << std::endl
		<< "z          = " << ph_space.z       << std::endl
		<< "ph_t²      = " << ph_space.ph_t_sq << std::endl
		<< "φ_h        = " << ph_space.phi_h   << std::endl
		<< "φ          = " << ph_space.phi     << std::endl
		<< "λ_e        = " << lambda         << std::endl
		<< "η_1        = " << eta.x            << std::endl
		<< "η_2        = " << eta.y            << std::endl
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

