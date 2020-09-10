#include <catch2/catch.hpp>

#include <fstream>
#include <istream>
#include <limits>
#include <sstream>
#include <utility>

#include <sidis/constant.hpp>
#include <sidis/kinematics.hpp>
#include <sidis/math.hpp>
#include <sidis/vector.hpp>

#include "stream_generator.hpp"
#include "rel_matcher.hpp"

using namespace sidis;

namespace {

struct Input {
	char particle_id;
	Real beam_energy;
	kin::PhaseSpace phase_space;
};

std::istream& operator>>(std::istream& in, Input& input) {
	in >> input.particle_id;
	if (!(input.particle_id == 'e'
			|| input.particle_id == 'm'
			|| input.particle_id == 't')) {
		in.setstate(std::ios_base::failbit);
	}
	in >> input.beam_energy;
	in >> input.phase_space.x;
	in >> input.phase_space.y;
	in >> input.phase_space.z;
	in >> input.phase_space.ph_t_sq;
	in >> input.phase_space.phi_h;
	in >> input.phase_space.phi;
	return in;
}

}

TEST_CASE(
		"Non-radiative kinematics checks",
		"[kin]") {
	Input input = GENERATE(
		from_stream<Input>(
			std::move(std::ifstream("data/phase_space_vals.dat")),
			true));
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
	kin::Initial initial_state(M, m, E_b);
	kin::PhaseSpace phase_space = input.phase_space;
	kin::Kinematics kin(initial_state, phase_space, mh, M_th);
	kin::Final final_state(initial_state, kin);
	// Get 4-momenta of particles.
	math::Vec4 p = initial_state.p;
	math::Vec4 k1 = initial_state.k1;
	math::Vec4 q = final_state.q;
	math::Vec4 k2 = final_state.k2;
	math::Vec4 ph = final_state.ph;
	math::Vec4 px = (p + k1) - (k2 + ph);
	// Basis vectors for angle checks.
	math::Vec3 e_y = cross(k1.r, k2.r).unit();
	math::Vec3 e_x = cross(e_y, q.r).unit();

	// Print state information.
	std::stringstream ss;
	ss
		<< "pid   = " << input.particle_id   << std::endl
		<< "E_b   = " << E_b                 << std::endl
		<< "x     = " << phase_space.x       << std::endl
		<< "y     = " << phase_space.y       << std::endl
		<< "z     = " << phase_space.z       << std::endl
		<< "ph_t² = " << phase_space.ph_t_sq << std::endl
		<< "φ_h   = " << phase_space.phi_h   << std::endl
		<< "φ     = " << phase_space.phi     << std::endl;
	INFO(ss.str());

	// Do comparisons.
	Real prec = 1e4 * std::numeric_limits<Real>::epsilon();

	// Kinematic variables.
	CHECK_THAT(
		-dot(q, q)/(2.*dot(q, p)),
		RelMatcher<Real>(kin.x, prec));
	CHECK_THAT(
		dot(q, p)/dot(k1, p),
		RelMatcher<Real>(kin.y, prec));
	CHECK_THAT(
		dot(ph, p)/dot(p, q),
		RelMatcher<Real>(kin.z, prec));
	CHECK_THAT(
		(q - ph).norm_sq(),
		RelMatcher<Real>(kin.t, prec));
	CHECK_THAT(
		2.*dot(p, k1),
		RelMatcher<Real>(kin.S, prec));
	CHECK_THAT(
		-q.norm_sq(),
		RelMatcher<Real>(kin.Q_sq, prec));
	CHECK_THAT(
		2.*dot(p, k2),
		RelMatcher<Real>(kin.X, prec));
	CHECK_THAT(
		2.*dot(p, q),
		RelMatcher<Real>(kin.S_x, prec));
	CHECK_THAT(
		2.*dot(k1, ph),
		RelMatcher<Real>(kin.V_1, prec));
	CHECK_THAT(
		2.*dot(k2, ph),
		RelMatcher<Real>(kin.V_2, prec));
	CHECK_THAT(
		dot(q, ph),
		RelMatcher<Real>(kin.V_m, prec));

	// 3-momenta magnitudes.
	CHECK_THAT(
		k1.r.norm(),
		RelMatcher<Real>(kin.lambda_S_sqrt/(2.*M), prec));
	CHECK_THAT(
		q.r.norm(),
		RelMatcher<Real>(kin.lambda_Y_sqrt/(2.*M), prec));

	// 3-momenta components.
	CHECK_THAT(
		dot(ph.r, q.r.unit()),
		RelMatcher<Real>(kin.ph_l, prec));
	CHECK_THAT(
		cross(ph.r, q.r.unit()).norm(),
		RelMatcher<Real>(kin.ph_t, prec));
	CHECK_THAT(
		cross(k1.r, q.r.unit()).norm(),
		RelMatcher<Real>(kin.k_t, prec));
	CHECK_THAT(
		cross(k2.r, q.r.unit()).norm(),
		RelMatcher<Real>(kin.k_t, prec));

	// Volume part.
	CHECK_THAT(
		dot(cross(k1.r, q.r), ph.r),
		RelMatcher<Real>(kin.vol_phi_h/M, prec));

	// Angles.
	CHECK_THAT(
		std::atan2(dot(e_y, ph.r), dot(e_x, ph.r)),
		RelMatcher<Real>(kin.phi_h, prec));
	CHECK_THAT(
		std::atan2(k2.r.y, k2.r.x),
		RelMatcher<Real>(kin.phi, prec));

	// Particle masses.
	Real prec_base = 1e4 * std::numeric_limits<Real>::epsilon();
	// The precisions for the comparisons are calculated in this way because
	// often there is a lot of precision lost when calculating the mass of a
	// particle with high energy.
	CHECK_THAT(
		p.norm(),
		RelMatcher<Real>(M, prec));
	CHECK_THAT(
		k1.norm(),
		RelMatcher<Real>(m, prec_base / (1. - k1.r.norm_sq()/math::sq(k1.t))));
	CHECK_THAT(
		k2.norm(),
		RelMatcher<Real>(m, prec_base / (1. - k2.r.norm_sq()/math::sq(k2.t))));
	CHECK_THAT(
		ph.norm(),
		RelMatcher<Real>(mh, prec_base / (1. - ph.r.norm_sq()/math::sq(ph.t))));
	CHECK_THAT(
		px.norm(),
		RelMatcher<Real>(kin.mx, prec_base / (1. - px.r.norm_sq()/math::sq(px.t))));
}

