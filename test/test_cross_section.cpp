#include <catch2/catch.hpp>

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
};

}

TEST_CASE(
		"Cross-section values",
		"[xs]") {
	// Load the structure function data ahead of time.
	sf::model::WW ww;

	using Tuple = std::tuple<Input, Output>;
	// Randomly generated cross-sections points for comparison.
	auto target = GENERATE(
		Tuple {
			{
				272.3289320834889,
				0.7742207157181373,
				{ 0.0100723080218603, -0.0074201560094463, -0.0861055915456175 },
				{ 0.2438003952188158, 0.1172111539548094, 0.7003402901345090, 0.9374573911481556, -0.5714773016683213, -2.179063668473446 },
			},
			{
				RelMatcher<Real>(4.95243870781e-7, 1e-11),
				RelMatcher<Real>(9.6695598884e-18, 1e-10),
				RelMatcher<Real>(-96.666414580, 1e-11),
				RelMatcher<Real>(14.8883577396684, 1e-14),
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
	// Choose beam and target polarizations.
	Real beam_pol = input.beam_pol;
	math::Vec3 target_pol = input.target_pol;
	// Transform the target polarization into a different basis.
	math::Vec3 e_z = final_state.q.r.unit();
	math::Vec3 e_y = math::cross(e_z, final_state.ph.r).unit();
	math::Vec3 e_x = math::cross(e_y, e_z).unit();
	math::Vec3 eta(
		math::dot(target_pol, e_x),
		math::dot(target_pol, e_y),
		math::dot(target_pol, e_z));
	// Calculate structure functions.
	sf::Sf sf = ww.sf(kin.x, kin.z, kin.Q_sq, kin.ph_t);
	// Compute the cross-sections.
	Real born = xs::born(beam_pol, eta, kin, sf);
	Real amm = xs::amm(beam_pol, eta, kin, sf);
	Real delta_vr = xs::delta_vr(kin);
	Real delta_vac_lep = xs::delta_vac_lep(kin);

	// Do comparisons.
	CHECK_THAT(born, output.born);
	CHECK_THAT(amm, output.amm);
	CHECK_THAT(delta_vr, output.delta_vr);
	CHECK_THAT(delta_vac_lep, output.delta_vac_lep);
}

TEST_CASE(
		"Cross-section values with muon",
		"[xs]") {
	// Test in the regime where `m` is comparable to `M`.
}

