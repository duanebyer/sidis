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
		},
		Tuple {
			{
				195.3294172569680,
				0.9100140348590521,
				{ -0.08802730180322475, -0.2533988547013961, 0.1187555736304225 },
				{ 0.4306540661961004, 0.1158543745257551, 0.9686469439861846, 0.02139603856591088, 2.192299006199159, -0.4261416174684685 },
			},
			{
				RelMatcher<Real>(3.5179877343e-9, 3.3e-11),
				RelMatcher<Real>(4.520130460e-20, 5.0e-10),
				RelMatcher<Real>(-195.716972083, 6.3e-12),
				RelMatcher<Real>(15.2623319875021, 3.3e-15),
				RelMatcher<Real>(8.94937158518684, 1.0e-15),
			},
		},
		Tuple {
			{
				62.59175892888381,
				-0.3275983794705739,
				{ -0.2047597459728121, -0.4040812880364820, -0.8626351840758933 },
				{ 0.4389898502347507, 0.4124597479558581, 0.5361865050334323, 0.7449474723556813, -0.1842149526004252, -0.4698043138784573 },
			},
			{
				RelMatcher<Real>(6.92450895176e-8, 1.6e-12),
				RelMatcher<Real>(1.734752835519e-17, 1.5e-12),
				RelMatcher<Real>(-45.033908839, 9.2e-12),
				RelMatcher<Real>(15.5177996600907, 3.3e-15),
				RelMatcher<Real>(9.30141401789073, 1.0e-15),
			},
		},
		Tuple {
			{
				91.64036306921523,
				-0.8489411406153405,
				{ -0.1594604964413398, -0.02044235379020882, 0.02103776788066295 },
				{ 0.2276654751955528, 0.5307882077424426, 0.5860674210347391, 0.8834491416986462, -1.480372340506740, 2.564777606223223 },
			},
			{
				RelMatcher<Real>(2.32005737013e-7, 3.8e-12),
				RelMatcher<Real>(5.58280288634e-17, 6.8e-13),
				RelMatcher<Real>(-21.071819341, 4.4e-11),
				RelMatcher<Real>(15.4783970324341, 3.3e-15),
				RelMatcher<Real>(9.24731032694575, 1.0e-15),
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

