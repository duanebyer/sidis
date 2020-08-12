#include <catch2/catch.hpp>

#include <sstream>
#include <tuple>

#include <sidis/sidis.hpp>
#include <sidis/sf_model/ww.hpp>

#include "rel_matcher.hpp"

using namespace sidis;

namespace {

struct Input {
	Real beam_energy;
	Real beam_pol;
	math::Vec3 target_pol;
	kin::PhaseSpace phase_space;
};

struct Output {
	RelMatcher<Real> born;
	RelMatcher<Real> amm;
	RelMatcher<Real> delta_vr;
	RelMatcher<Real> delta_vac_lep;
	RelMatcher<Real> delta_vac_had;
};

}

TEST_CASE(
		"Non-radiative cross-section values",
		"[xs]") {
	// Load the structure function data ahead of time.
	sf::model::WW ww;

	using Tuple = std::tuple<Input, Output>;
	// Randomly generated cross-sections points for comparison.
	auto target = GENERATE(
		Tuple {
			{
				244.3994717676956,
				-0.6164827543415216,
				{ -0.4557051628146152, 0.09688945772348936, -0.4280942449984859 },
				{ 0.2195446560618665, 0.03029886835408146, 0.2387660175848276, 0.5575045307531898, -2.344259584082282, 0.1116601857460971 },
			},
			{
				RelMatcher<Real>(1.2860080823885e-3, 2.1e-13),
				RelMatcher<Real>(1.97833608719e-14, 3.1e-12),
				RelMatcher<Real>(-99.4643915906, 1.0e-12),
				RelMatcher<Real>(12.4986683926813, 4.2e-15),
				RelMatcher<Real>(5.07031143940847, 1.2e-15),
			},
		});
	Input input = std::get<0>(target);
	Output output = std::get<1>(target);

	// Set up the input to the cross-section calculation.
	Real E_b = input.beam_energy;
	Real M = constant::MASS_P;
	Real m = constant::MASS_E;
	Real mh = constant::MASS_PI;
	Real M_th = constant::MASS_P + constant::MASS_PI_0;
	Real pi = constant::PI;
	kin::Initial initial_state {
		math::Vec4(M, 0., 0., 0.),
		math::Vec4(E_b, std::sqrt(E_b * E_b - m * m) * math::Vec3::Z),
	};
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
	ss << "x     = " << std::to_string(phase_space.x)       << std::endl;
	ss << "y     = " << std::to_string(phase_space.y)       << std::endl;
	ss << "z     = " << std::to_string(phase_space.z)       << std::endl;
	ss << "ph_t² = " << std::to_string(phase_space.ph_t_sq) << std::endl;
	ss << "φ_h   = " << std::to_string(phase_space.phi_h)   << std::endl;
	ss << "φ     = " << std::to_string(phase_space.phi)     << std::endl;
	ss << "λ_e   = " << std::to_string(beam_pol)            << std::endl;
	ss << "η_1   = " << std::to_string(target_pol.x)        << std::endl;
	ss << "η_2   = " << std::to_string(target_pol.y)        << std::endl;
	ss << "η_3   = " << std::to_string(target_pol.z);
	INFO(ss.str());

	// Do comparisons.
	CHECK_THAT(born, output.born);
	CHECK_THAT(amm, output.amm);
	CHECK_THAT(delta_vr, output.delta_vr);
	CHECK_THAT(delta_vac_lep, output.delta_vac_lep);
	CHECK_THAT(delta_vac_had, output.delta_vac_had);
}

TEST_CASE(
		"Non-radiative cross-section values with muon",
		"[xs]") {
	// Test in the regime where `m` is comparable to `M`.
}

