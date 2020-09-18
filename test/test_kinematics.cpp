#include <catch2/catch.hpp>

#include <fstream>
#include <istream>
#include <limits>
#include <sstream>
#include <utility>

#include <sidis/constant.hpp>
#include <sidis/kinematics.hpp>
#include <sidis/extra/math.hpp>
#include <sidis/extra/vector.hpp>

#include "rel_matcher.hpp"
#include "stream_generator.hpp"

using namespace sidis;

namespace {

struct Input {
	char particle_id;
	Real beam_energy;
	kin::PhaseSpace phase_space;
};

struct InputRad {
	char particle_id;
	Real beam_energy;
	kin::PhaseSpaceRad phase_space;
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

std::istream& operator>>(std::istream& in, InputRad& input) {
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
	in >> input.phase_space.tau;
	in >> input.phase_space.R;
	in >> input.phase_space.phi_k;
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
	math::Vec3 e_y = cross(k1.vec3(), k2.vec3()).unit();
	math::Vec3 e_x = cross(e_y, q.vec3()).unit();

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
		k1.vec3().norm(),
		RelMatcher<Real>(kin.lambda_S_sqrt/(2.*M), prec));
	CHECK_THAT(
		q.vec3().norm(),
		RelMatcher<Real>(kin.lambda_Y_sqrt/(2.*M), prec));

	// 3-momenta components.
	CHECK_THAT(
		dot(ph.vec3(), q.vec3().unit()),
		RelMatcher<Real>(kin.ph_l, prec));
	CHECK_THAT(
		cross(ph.vec3(), q.vec3().unit()).norm(),
		RelMatcher<Real>(kin.ph_t, prec));
	CHECK_THAT(
		cross(k1.vec3(), q.vec3().unit()).norm(),
		RelMatcher<Real>(kin.k1_t, prec));
	CHECK_THAT(
		cross(k2.vec3(), q.vec3().unit()).norm(),
		RelMatcher<Real>(kin.k1_t, prec));

	// Volume parts.
	CHECK_THAT(
		dot(cross(k1.vec3(), q.vec3()), ph.vec3()),
		RelMatcher<Real>(kin.vol_phi_h/M, prec));

	// Angles.
	CHECK_THAT(
		std::atan2(dot(e_y, ph.vec3()), dot(e_x, ph.vec3())),
		RelMatcher<Real>(kin.phi_h, prec));
	CHECK_THAT(
		std::atan2(k2.vec3().y, k2.vec3().x),
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
		RelMatcher<Real>(m, prec_base / (1. - k1.vec3().norm_sq()/math::sq(k1.t))));
	CHECK_THAT(
		k2.norm(),
		RelMatcher<Real>(m, prec_base / (1. - k2.vec3().norm_sq()/math::sq(k2.t))));
	CHECK_THAT(
		ph.norm(),
		RelMatcher<Real>(mh, prec_base / (1. - ph.vec3().norm_sq()/math::sq(ph.t))));
	CHECK_THAT(
		px.norm(),
		RelMatcher<Real>(kin.mx, prec_base / (1. - px.vec3().norm_sq()/math::sq(px.t))));
}

TEST_CASE(
		"Radiative kinematics checks",
		"[kin]") {
	InputRad input = GENERATE(
		from_stream<InputRad>(
			std::move(std::ifstream("data/phase_space_rad_vals.dat")),
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
	kin::PhaseSpaceRad phase_space = input.phase_space;
	kin::KinematicsRad kin(initial_state, phase_space, mh, M_th);
	kin::FinalRad final_state(initial_state, kin);
	// Get 4-momenta of particles.
	math::Vec4 p = initial_state.p;
	math::Vec4 k1 = initial_state.k1;
	math::Vec4 q = final_state.q;
	math::Vec4 k2 = final_state.k2;
	math::Vec4 k = final_state.k;
	math::Vec4 ph = final_state.ph;
	math::Vec4 px = (p + k1) - (k2 + ph);
	// Basis vectors for angle checks.
	math::Vec3 e_y = cross(k1.vec3(), k2.vec3()).unit();
	math::Vec3 e_x = cross(e_y, q.vec3()).unit();
	math::Vec3 shift_e_y = cross(k1.vec3(), (k2 + k).vec3()).unit();
	math::Vec3 shift_e_x = cross(shift_e_y, (q - k).vec3()).unit();

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
		<< "φ     = " << phase_space.phi     << std::endl
		<< "τ     = " << phase_space.tau     << std::endl
		<< "R     = " << phase_space.R       << std::endl
		<< "φ_k   = " << phase_space.phi_k   << std::endl;
	INFO(ss.str());

	// Do comparisons.
	Real prec = 1e4 * std::numeric_limits<Real>::epsilon();

	// Kinematic variables.
	CHECK_THAT(
		dot(k, q) / dot(k, p),
		RelMatcher<Real>(kin.tau, prec));
	CHECK_THAT(
		2. * dot(k, p),
		RelMatcher<Real>(kin.R, prec));
	CHECK_THAT(
		dot(k1, k) / dot(p, k),
		RelMatcher<Real>(kin.z_1, prec));
	CHECK_THAT(
		dot(k2, k) / dot(p, k),
		RelMatcher<Real>(kin.z_2, prec));
	CHECK_THAT(
		dot(k, ph) / dot(k, p),
		RelMatcher<Real>(kin.mu, prec));

	// 3-momenta dot products.
	CHECK_THAT(
		dot(k.vec3(), ph.vec3()),
		RelMatcher<Real>(kin.lambda_RV/(4.*M*M), prec));
	CHECK_THAT(
		dot(k.vec3(), q.vec3()),
		RelMatcher<Real>(kin.lambda_RY/(4.*M*M), prec));
	CHECK_THAT(
		dot(q.vec3(), ph.vec3()),
		RelMatcher<Real>(kin.lambda_V/(4.*M*M), prec));
	CHECK_THAT(
		ph.vec3().norm_sq(),
		RelMatcher<Real>(kin.lambda_H/(4.*M*M), prec));

	// Volume parts.
	CHECK_THAT(
		dot(cross(k1.vec3(), q.vec3()), k.vec3()),
		RelMatcher<Real>(kin.vol_phi_k/M, prec));
	CHECK_THAT(
		dot(cross(k.vec3(), ph.vec3()), q.vec3()),
		RelMatcher<Real>(kin.vol_phi_hk/M, prec));

	// Angles.
	CHECK_THAT(
		std::atan2(dot(e_y, k.vec3()), dot(e_x, k.vec3())),
		RelMatcher<Real>(kin.phi_k, prec));

	// Particle masses.
	CHECK(
		std::abs(k.norm_sq()) <=
		std::abs(prec * k.t * k.t));

	// Shifted kinematic variables.
	CHECK_THAT(
		-dot(q - k, q - k)/(2.*dot(q - k, p)),
		RelMatcher<Real>(kin.shift_x, prec));
	CHECK_THAT(
		dot(ph, p)/dot(p, q - k),
		RelMatcher<Real>(kin.shift_z, prec));
	CHECK_THAT(
		(q - k - ph).norm_sq(),
		RelMatcher<Real>(kin.shift_t, prec));
	CHECK_THAT(
		-(q - k).norm_sq(),
		RelMatcher<Real>(kin.shift_Q_sq, prec));
	CHECK_THAT(
		2.*dot(p, q - k),
		RelMatcher<Real>(kin.shift_S_x, prec));
	CHECK_THAT(
		dot(q - k, ph),
		RelMatcher<Real>(kin.shift_V_m, prec));

	// Shifted 3-momenta magnitudes.
	CHECK_THAT(
		(q - k).vec3().norm(),
		RelMatcher<Real>(kin.shift_lambda_Y_sqrt/(2.*M), prec));

	// Shifted 3-momenta components.
	CHECK_THAT(
		dot(ph.vec3(), (q - k).vec3().unit()),
		RelMatcher<Real>(kin.shift_ph_l, prec));
	CHECK_THAT(
		cross(ph.vec3(), (q - k).vec3().unit()).norm(),
		RelMatcher<Real>(kin.shift_ph_t, prec));
	CHECK_THAT(
		cross(k1.vec3(), (q - k).vec3().unit()).norm(),
		RelMatcher<Real>(kin.shift_k1_t, prec));

	// Shifted volume parts.
	CHECK_THAT(
		dot(cross(k1.vec3(), (q - k).vec3()), ph.vec3()),
		RelMatcher<Real>(kin.shift_vol_phi_h/M, prec));

	// Shifted angles.
	CHECK_THAT(
		std::atan2(dot(shift_e_y, ph.vec3()), dot(shift_e_x, ph.vec3())),
		RelMatcher<Real>(kin.shift_phi_h, prec));
	CHECK_THAT(
		std::atan2(k2.vec3().y, k2.vec3().x),
		RelMatcher<Real>(kin.phi, prec));
}

