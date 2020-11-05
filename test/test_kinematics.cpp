#include <catch2/catch.hpp>

#include <fstream>
#include <istream>
#include <limits>
#include <sstream>
#include <utility>

#include <sidis/constant.hpp>
#include <sidis/frame.hpp>
#include <sidis/kinematics.hpp>
#include <sidis/extra/math.hpp>
#include <sidis/extra/transform.hpp>
#include <sidis/extra/vector.hpp>

#include "abs_matcher.hpp"
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

void test_kin_nrad(
	kin::Initial initial_state,
	kin::Kinematics kin,
	bool complete);

TEST_CASE(
		"Non-radiative kinematics checks",
		"[kin]") {
	Input input = GENERATE(
		from_stream<Input>(
			std::move(std::ifstream("data/phase_space_vals.dat")),
			true));
	Real E_b = input.beam_energy;
	Real M = constant::MASS_P;
	constant::Lepton lep;
	if (input.particle_id == 'e') {
		lep = constant::Lepton::E;
	} else if (input.particle_id == 'm') {
		lep = constant::Lepton::MU;
	} else if (input.particle_id == 't') {
		lep = constant::Lepton::TAU;
	}
	Real M_th = constant::MASS_P + constant::MASS_PI_0;
	kin::Initial initial_state(constant::Nucleus::P, lep, E_b);
	kin::PhaseSpace phase_space = input.phase_space;
	kin::Kinematics kin(initial_state, phase_space, constant::Hadron::PI_P, M_th);

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

	test_kin_nrad(initial_state, kin, true);
}

TEST_CASE(
		"Non-radiative shifted kinematics checks",
		"[kin]") {
	InputRad input = GENERATE(
		from_stream<InputRad>(
			std::move(std::ifstream("data/phase_space_rad_vals.dat")),
			true));
	Real E_b = input.beam_energy;
	Real M = constant::MASS_P;
	constant::Lepton lep;
	if (input.particle_id == 'e') {
		lep = constant::Lepton::E;
	} else if (input.particle_id == 'm') {
		lep = constant::Lepton::MU;
	} else if (input.particle_id == 't') {
		lep = constant::Lepton::TAU;
	}
	Real M_th = constant::MASS_P + constant::MASS_PI_0;
	kin::Initial initial_state(constant::Nucleus::P, lep, E_b);
	kin::PhaseSpaceRad phase_space = input.phase_space;
	kin::KinematicsRad kin(initial_state, phase_space, constant::Hadron::PI_P, M_th);

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

	test_kin_nrad(initial_state, kin.project_shift(), false);
}

void test_kin_nrad(
		kin::Initial initial_state,
		kin::Kinematics kin,
		bool complete) {
	kin::Final final_state(initial_state, math::Vec3::Y, kin);
	// Get 4-momenta of particles.
	math::Vec4 p = initial_state.p;
	math::Vec4 k1 = initial_state.k1;
	math::Vec4 q = final_state.q;
	math::Vec4 k2 = final_state.k2;
	math::Vec4 ph = final_state.ph;
	math::Vec4 px = (p + k1) - (k2 + ph);
	// Basis vectors for angle checks.
	math::Vec3 e_y = cross(q.r(), k1.r()).unit();
	math::Vec3 e_x = cross(e_y, q.r()).unit();

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
		k1.r().norm(),
		RelMatcher<Real>(kin.lambda_S_sqrt/(2.*kin.M), prec));
	CHECK_THAT(
		q.r().norm(),
		RelMatcher<Real>(kin.lambda_Y_sqrt/(2.*kin.M), prec));

	// 3-momenta components.
	CHECK_THAT(
		dot(ph.r(), q.r().unit()),
		RelMatcher<Real>(kin.ph_l, prec));
	CHECK_THAT(
		cross(ph.r(), q.r().unit()).norm(),
		RelMatcher<Real>(kin.ph_t, prec));
	CHECK_THAT(
		cross(k1.r(), q.r().unit()).norm(),
		RelMatcher<Real>(kin.k1_t, prec));
	if (complete) {
		CHECK_THAT(
			cross(k2.r(), q.r().unit()).norm(),
			RelMatcher<Real>(kin.k1_t, prec));
	}
	CHECK_THAT(
		dot(q.r(), k1.r().unit()),
		RelMatcher<Real>(kin.q_l, prec));
	CHECK_THAT(
		cross(q.r(), k1.r().unit()).norm(),
		RelMatcher<Real>(kin.q_t, prec));

	// Volume parts.
	CHECK_THAT(
		dot(cross(k1.r(), q.r()), ph.r()),
		RelMatcher<Real>(kin.vol_phi_h/kin.M, prec));

	// Angles.
	CHECK_THAT(
		std::atan2(dot(e_y, ph.r()), dot(e_x, ph.r())),
		RelMatcher<Real>(kin.phi_h, prec));
	CHECK_THAT(
		std::atan2(-k2.x, k2.y),
		RelMatcher<Real>(kin.phi, prec));
	CHECK_THAT(
		std::atan2(-q.x, q.y),
		RelMatcher<Real>(kin.phi_q, prec));

	// Completeness.
	if (complete) {
		math::Vec4 q_test = k1 - k2;
		CHECK_THAT(
			q_test.t,
			RelMatcher<Real>(q.t, prec));
		CHECK_THAT(
			q_test.x,
			RelMatcher<Real>(q.x, prec));
		CHECK_THAT(
			q_test.y,
			RelMatcher<Real>(q.y, prec));
		CHECK_THAT(
			q_test.z,
			RelMatcher<Real>(q.z, prec));
	}

	// Particle masses.
	Real prec_base = 1e4 * std::numeric_limits<Real>::epsilon();
	// The precisions for the comparisons are calculated in this way because
	// often there is a lot of precision lost when calculating the mass of a
	// particle with high energy.
	CHECK_THAT(
		p.norm(),
		RelMatcher<Real>(kin.M, prec));
	CHECK_THAT(
		k1.norm(),
		RelMatcher<Real>(kin.m, prec_base/(1. - k1.r().norm_sq()/math::sq(k1.t))));
	CHECK_THAT(
		k2.norm(),
		RelMatcher<Real>(kin.m, prec_base/(1. - k2.r().norm_sq()/math::sq(k2.t))));
	CHECK_THAT(
		ph.norm(),
		RelMatcher<Real>(kin.mh, prec_base/(1. - ph.r().norm_sq()/math::sq(ph.t))));
	if (complete) {
		CHECK_THAT(
			px.norm(),
			RelMatcher<Real>(kin.mx, prec_base/(1. - px.r().norm_sq()/math::sq(px.t))));
	}
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
	constant::Lepton lep;
	if (input.particle_id == 'e') {
		lep = constant::Lepton::E;
	} else if (input.particle_id == 'm') {
		lep = constant::Lepton::MU;
	} else if (input.particle_id == 't') {
		lep = constant::Lepton::TAU;
	}
	Real M_th = constant::MASS_P + constant::MASS_PI_0;
	kin::Initial initial_state(constant::Nucleus::P, lep, E_b);
	kin::PhaseSpaceRad phase_space = input.phase_space;
	kin::KinematicsRad kin(initial_state, phase_space, constant::Hadron::PI_P, M_th);
	kin::FinalRad final_state(initial_state, math::Vec3::ZERO, kin);
	// Get 4-momenta of particles.
	math::Vec4 p = initial_state.p;
	math::Vec4 k1 = initial_state.k1;
	math::Vec4 q = final_state.q;
	math::Vec4 k2 = final_state.k2;
	math::Vec4 k = final_state.k;
	math::Vec4 ph = final_state.ph;
	math::Vec4 px = (p + k1) - (k2 + ph);
	// Basis vectors for angle checks.
	math::Vec3 e_y = cross(k1.r(), k2.r()).unit();
	math::Vec3 e_x = cross(e_y, q.r()).unit();
	math::Vec3 shift_e_y = cross(k1.r(), (k2 + k).r()).unit();
	math::Vec3 shift_e_x = cross(shift_e_y, (q - k).r()).unit();

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
		dot(k, q)/dot(k, p),
		RelMatcher<Real>(kin.tau, prec));
	CHECK_THAT(
		2.*dot(k, p),
		RelMatcher<Real>(kin.R, prec));
	CHECK_THAT(
		dot(k1, k)/dot(p, k),
		RelMatcher<Real>(kin.z_1, prec));
	CHECK_THAT(
		dot(k2, k)/dot(p, k),
		RelMatcher<Real>(kin.z_2, prec));
	CHECK_THAT(
		dot(k, ph)/dot(k, p),
		RelMatcher<Real>(kin.mu, prec));
	CHECK_THAT(
		(k1 / kin.z_1 - k2 / kin.z_2).norm_sq(),
		RelMatcher<Real>(kin.F_IR, prec));

	// 3-momenta dot products.
	CHECK_THAT(
		dot(k.r(), ph.r()),
		RelMatcher<Real>(kin.lambda_RV/(4.*M*M), prec));
	CHECK_THAT(
		dot(k.r(), q.r()),
		RelMatcher<Real>(kin.lambda_RY/(4.*M*M), prec));
	CHECK_THAT(
		dot(q.r(), ph.r()),
		RelMatcher<Real>(kin.lambda_V/(4.*M*M), prec));
	CHECK_THAT(
		ph.r().norm_sq(),
		RelMatcher<Real>(kin.lambda_H/(4.*M*M), prec));

	// Volume parts.
	CHECK_THAT(
		dot(cross(k1.r(), q.r()), k.r()),
		RelMatcher<Real>(kin.R*kin.vol_phi_k_R/M, prec));
	CHECK_THAT(
		dot(cross(k.r(), ph.r()), q.r()),
		RelMatcher<Real>(kin.vol_phi_hk/M, prec));

	// Angles.
	CHECK_THAT(
		std::atan2(dot(e_y, k.r()), dot(e_x, k.r())),
		RelMatcher<Real>(kin.phi_k, prec));

	// Particle masses.
	CHECK(
		std::abs(k.norm_sq()) <=
		std::abs(prec*k.t*k.t));

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
		(q - k).r().norm(),
		RelMatcher<Real>(kin.shift_lambda_Y_sqrt/(2.*M), prec));

	// Shifted 3-momenta components.
	CHECK_THAT(
		dot(ph.r(), (q - k).r().unit()),
		RelMatcher<Real>(kin.shift_ph_l, prec));
	CHECK_THAT(
		cross(ph.r(), (q - k).r().unit()).norm(),
		RelMatcher<Real>(kin.shift_ph_t, prec));
	CHECK_THAT(
		cross(k1.r(), (q - k).r().unit()).norm(),
		RelMatcher<Real>(kin.shift_k1_t, prec));
	CHECK_THAT(
		dot((q - k).r(), k1.r().unit()),
		RelMatcher<Real>(kin.shift_q_l, prec));
	CHECK_THAT(
		cross((q - k).r(), k1.r().unit()).norm(),
		RelMatcher<Real>(kin.shift_q_t, prec));

	// Shifted volume parts.
	CHECK_THAT(
		dot(cross(k1.r(), (q - k).r()), ph.r()),
		RelMatcher<Real>(kin.shift_vol_phi_h/M, prec));

	// Shifted angles.
	CHECK_THAT(
		std::atan2(dot(shift_e_y, ph.r()), dot(shift_e_x, ph.r())),
		RelMatcher<Real>(kin.shift_phi_h, prec));
	CHECK_THAT(
		std::atan2((q - k).x, -(q - k).y),
		RelMatcher<Real>(kin.shift_phi_q, prec));
}

TEST_CASE(
		"Target to lab frame checks",
		"[frame]") {
	math::Vec3 p(1.2, -0.5, 2.3);
	math::Vec3 k1(0.4, -2.2, 0.2);
	math::Vec3 pol(0.2, 0.1, -0.3);
	kin::Initial initial_state(constant::Nucleus::P, p, constant::Lepton::MU, k1);

	// Construct the frames.
	math::Transform4 target_from_lab = frame::target_from_lab(initial_state, pol);
	math::Transform4 lab_from_target = frame::lab_from_target(initial_state, pol);
	math::Transform4 lab_from_lab = lab_from_target * target_from_lab;

	Real prec = 1e4 * std::numeric_limits<Real>::epsilon();

	// Check orthogonality.
	CHECK_THAT(lab_from_lab.t.t, AbsMatcher<Real>(1., prec));
	CHECK_THAT(lab_from_lab.t.x, AbsMatcher<Real>(0., prec));
	CHECK_THAT(lab_from_lab.t.y, AbsMatcher<Real>(0., prec));
	CHECK_THAT(lab_from_lab.t.z, AbsMatcher<Real>(0., prec));
	CHECK_THAT(lab_from_lab.x.t, AbsMatcher<Real>(0., prec));
	CHECK_THAT(lab_from_lab.x.x, AbsMatcher<Real>(-1., prec));
	CHECK_THAT(lab_from_lab.x.y, AbsMatcher<Real>(0., prec));
	CHECK_THAT(lab_from_lab.x.z, AbsMatcher<Real>(0., prec));
	CHECK_THAT(lab_from_lab.y.t, AbsMatcher<Real>(0., prec));
	CHECK_THAT(lab_from_lab.y.x, AbsMatcher<Real>(0., prec));
	CHECK_THAT(lab_from_lab.y.y, AbsMatcher<Real>(-1., prec));
	CHECK_THAT(lab_from_lab.y.z, AbsMatcher<Real>(0., prec));
	CHECK_THAT(lab_from_lab.z.t, AbsMatcher<Real>(0., prec));
	CHECK_THAT(lab_from_lab.z.x, AbsMatcher<Real>(0., prec));
	CHECK_THAT(lab_from_lab.z.y, AbsMatcher<Real>(0., prec));
	CHECK_THAT(lab_from_lab.z.z, AbsMatcher<Real>(-1., prec));

	// Check determinants.
	CHECK_THAT(target_from_lab.det(), RelMatcher<Real>(1., prec));
	CHECK_THAT(lab_from_target.det(), RelMatcher<Real>(1., prec));

	// Check transforming `p` and `k1` correctly.
	math::Vec4 p_target = target_from_lab * initial_state.p;
	math::Vec4 k1_target = target_from_lab * initial_state.k1;
	CHECK_THAT(p_target.t, RelMatcher<Real>(constant::MASS_P, prec));
	CHECK_THAT(p_target.x, AbsMatcher<Real>(0., prec));
	CHECK_THAT(p_target.y, AbsMatcher<Real>(0., prec));
	CHECK_THAT(p_target.z, AbsMatcher<Real>(0., prec));
	CHECK(k1_target.t > 0.);
	CHECK_THAT(k1_target.x, AbsMatcher<Real>(0., prec));
	CHECK_THAT(k1_target.y, AbsMatcher<Real>(0., prec));
	CHECK(k1_target.z > 0.);
	CHECK_THAT(k1_target.norm(), RelMatcher<Real>(constant::MASS_MU, prec));

	// Check transforming polarization correctly. Remember that polarization is
	// given in the proton rest frame, so first have to boost it into proton
	// frame.
	math::Transform4 pol_bv(
		0., 0., 0., 0.,
		0., 0., -pol.z, pol.y,
		0., pol.z, 0., -pol.x,
		0., -pol.y, pol.x, 0.);
	math::Transform4 boost = math::Transform4::boost_to(initial_state.p);
	math::Vec4 k1_boost = boost.transpose() * initial_state.k1;
	math::Transform4 pol_lab_bv = boost.transform(pol_bv);
	math::Transform4 pol_target_bv = target_from_lab.transform(pol_lab_bv);
	math::Vec3 pol_target(
		pol_target_bv.y.z,
		pol_target_bv.z.x,
		pol_target_bv.x.y);
	Real pol_perp = pol.perp(k1_boost.r()).norm();
	Real pol_par = dot(pol, k1_boost.r().unit());
	CHECK_THAT(pol_target_bv.t.t, AbsMatcher<Real>(0., prec));
	CHECK_THAT(pol_target_bv.t.x, AbsMatcher<Real>(0., prec));
	CHECK_THAT(pol_target_bv.t.y, AbsMatcher<Real>(0., prec));
	CHECK_THAT(pol_target_bv.t.z, AbsMatcher<Real>(0., prec));
	CHECK_THAT(pol_target_bv.x.t, AbsMatcher<Real>(0., prec));
	CHECK_THAT(pol_target_bv.x.x, AbsMatcher<Real>(0., prec));
	CHECK_THAT(pol_target_bv.x.y, RelMatcher<Real>(pol_par, prec));
	CHECK_THAT(pol_target_bv.x.z, RelMatcher<Real>(-pol_perp, prec));
	CHECK_THAT(pol_target_bv.y.t, AbsMatcher<Real>(0., prec));
	CHECK_THAT(pol_target_bv.y.x, RelMatcher<Real>(-pol_par, prec));
	CHECK_THAT(pol_target_bv.y.y, AbsMatcher<Real>(0., prec));
	CHECK_THAT(pol_target_bv.y.z, AbsMatcher<Real>(0., prec));
	CHECK_THAT(pol_target_bv.z.t, AbsMatcher<Real>(0., prec));
	CHECK_THAT(pol_target_bv.z.x, RelMatcher<Real>(pol_perp, prec));
	CHECK_THAT(pol_target_bv.z.y, AbsMatcher<Real>(0., prec));
	CHECK_THAT(pol_target_bv.z.z, AbsMatcher<Real>(0., prec));
	// Check that angle is preserved.
	CHECK_THAT(k1_boost.r().norm(), RelMatcher<Real>(k1_target.r().norm(), prec));
	CHECK_THAT(pol.norm(), RelMatcher<Real>(pol_target.norm(), prec));
	CHECK_THAT(
		dot(pol_target, k1_target.r()),
		RelMatcher<Real>(dot(pol, k1_boost.r()), prec));
}

TEST_CASE(
		"Default target to lab frame checks",
		"[frame]") {
	kin::Initial initial_state(constant::Nucleus::P, constant::Lepton::MU, 8.2);

	// Construct the frames.
	math::Vec3 pol = GENERATE(
		math::Vec3(0., 0., 0.),
		math::Vec3(0., 0., 1e-8),
		math::Vec3(0., 1e-9, 1e-8),
		math::Vec3(0.1, -0.2, 0.3),
		math::Vec3(0., 0., -0.4));
	math::Transform4 target_from_lab = frame::target_from_lab(initial_state, pol);
	math::Transform4 lab_from_target = frame::lab_from_target(initial_state, pol);
	math::Transform4 lab_from_lab = lab_from_target * target_from_lab;

	Real prec = 1e4 * std::numeric_limits<Real>::epsilon();

	// Check that the frames are identity except for a rotation around z-axis.
	CHECK_THAT(target_from_lab.t.t, AbsMatcher<Real>(1., prec));
	CHECK_THAT(target_from_lab.t.x, AbsMatcher<Real>(0., prec));
	CHECK_THAT(target_from_lab.t.y, AbsMatcher<Real>(0., prec));
	CHECK_THAT(target_from_lab.t.z, AbsMatcher<Real>(0., prec));
	CHECK_THAT(target_from_lab.x.t, AbsMatcher<Real>(0., prec));
	CHECK_THAT(target_from_lab.x.z, AbsMatcher<Real>(0., prec));
	CHECK_THAT(target_from_lab.y.t, AbsMatcher<Real>(0., prec));
	CHECK_THAT(target_from_lab.y.z, AbsMatcher<Real>(0., prec));
	CHECK_THAT(target_from_lab.z.t, AbsMatcher<Real>(0., prec));
	CHECK_THAT(target_from_lab.z.x, AbsMatcher<Real>(0., prec));
	CHECK_THAT(target_from_lab.z.y, AbsMatcher<Real>(0., prec));
	CHECK_THAT(target_from_lab.z.z, AbsMatcher<Real>(-1., prec));

	// Check determinants.
	CHECK_THAT(target_from_lab.det(), RelMatcher<Real>(1., prec));
	CHECK_THAT(lab_from_target.det(), RelMatcher<Real>(1., prec));
}

TEST_CASE(
		"Non-radiative reference frame checks",
		"[frame]") {
	Input input = GENERATE(
		from_stream<Input>(
			std::move(std::ifstream("data/phase_space_vals.dat")),
			true));

	Real E_b = input.beam_energy;
	Real M = constant::MASS_P;
	constant::Lepton lep;
	if (input.particle_id == 'e') {
		lep = constant::Lepton::E;
	} else if (input.particle_id == 'm') {
		lep = constant::Lepton::MU;
	} else if (input.particle_id == 't') {
		lep = constant::Lepton::TAU;
	}
	Real M_th = constant::MASS_P + constant::MASS_PI_0;
	kin::Initial initial_state(constant::Nucleus::P, lep, E_b);
	kin::PhaseSpace phase_space = input.phase_space;
	kin::Kinematics kin(initial_state, phase_space, constant::Hadron::PI_P, M_th);
	kin::Final final_state(initial_state, math::Vec3::Y, kin);
	// Reference frames.
	math::Transform4 target_from_lepton = frame::target_from_lepton(kin);
	math::Transform4 target_from_hadron = frame::target_from_hadron(kin);
	math::Transform4 target_from_virt_photon = frame::target_from_virt_photon(kin);
	// Get 4-momenta of particles.
	math::Vec4 p = initial_state.p;
	math::Vec4 k1 = initial_state.k1;
	math::Vec4 q = final_state.q;
	math::Vec4 k2 = final_state.k2;
	math::Vec4 ph = final_state.ph;

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

	Real prec = 1e4 * std::numeric_limits<Real>::epsilon();

	// Check orthogonality.
	math::Transform4 lepton_from_lepton =
		target_from_lepton.transpose() * target_from_lepton;
	math::Transform4 hadron_from_hadron =
		target_from_hadron.transpose() * target_from_hadron;
	math::Transform4 virt_photon_from_virt_photon =
		target_from_virt_photon.transpose() * target_from_virt_photon;

	CHECK_THAT(lepton_from_lepton.t.t, AbsMatcher<Real>(1., prec));
	CHECK_THAT(lepton_from_lepton.t.x, AbsMatcher<Real>(0., prec));
	CHECK_THAT(lepton_from_lepton.t.y, AbsMatcher<Real>(0., prec));
	CHECK_THAT(lepton_from_lepton.t.z, AbsMatcher<Real>(0., prec));
	CHECK_THAT(lepton_from_lepton.x.t, AbsMatcher<Real>(0., prec));
	CHECK_THAT(lepton_from_lepton.x.x, AbsMatcher<Real>(-1., prec));
	CHECK_THAT(lepton_from_lepton.x.y, AbsMatcher<Real>(0., prec));
	CHECK_THAT(lepton_from_lepton.x.z, AbsMatcher<Real>(0., prec));
	CHECK_THAT(lepton_from_lepton.y.t, AbsMatcher<Real>(0., prec));
	CHECK_THAT(lepton_from_lepton.y.x, AbsMatcher<Real>(0., prec));
	CHECK_THAT(lepton_from_lepton.y.y, AbsMatcher<Real>(-1., prec));
	CHECK_THAT(lepton_from_lepton.y.z, AbsMatcher<Real>(0., prec));
	CHECK_THAT(lepton_from_lepton.z.t, AbsMatcher<Real>(0., prec));
	CHECK_THAT(lepton_from_lepton.z.x, AbsMatcher<Real>(0., prec));
	CHECK_THAT(lepton_from_lepton.z.y, AbsMatcher<Real>(0., prec));
	CHECK_THAT(lepton_from_lepton.z.z, AbsMatcher<Real>(-1., prec));

	CHECK_THAT(hadron_from_hadron.t.t, AbsMatcher<Real>(1., prec));
	CHECK_THAT(hadron_from_hadron.t.x, AbsMatcher<Real>(0., prec));
	CHECK_THAT(hadron_from_hadron.t.y, AbsMatcher<Real>(0., prec));
	CHECK_THAT(hadron_from_hadron.t.z, AbsMatcher<Real>(0., prec));
	CHECK_THAT(hadron_from_hadron.x.t, AbsMatcher<Real>(0., prec));
	CHECK_THAT(hadron_from_hadron.x.x, AbsMatcher<Real>(-1., prec));
	CHECK_THAT(hadron_from_hadron.x.y, AbsMatcher<Real>(0., prec));
	CHECK_THAT(hadron_from_hadron.x.z, AbsMatcher<Real>(0., prec));
	CHECK_THAT(hadron_from_hadron.y.t, AbsMatcher<Real>(0., prec));
	CHECK_THAT(hadron_from_hadron.y.x, AbsMatcher<Real>(0., prec));
	CHECK_THAT(hadron_from_hadron.y.y, AbsMatcher<Real>(-1., prec));
	CHECK_THAT(hadron_from_hadron.y.z, AbsMatcher<Real>(0., prec));
	CHECK_THAT(hadron_from_hadron.z.t, AbsMatcher<Real>(0., prec));
	CHECK_THAT(hadron_from_hadron.z.x, AbsMatcher<Real>(0., prec));
	CHECK_THAT(hadron_from_hadron.z.y, AbsMatcher<Real>(0., prec));
	CHECK_THAT(hadron_from_hadron.z.z, AbsMatcher<Real>(-1., prec));

	CHECK_THAT(virt_photon_from_virt_photon.t.t, AbsMatcher<Real>(1., prec));
	CHECK_THAT(virt_photon_from_virt_photon.t.x, AbsMatcher<Real>(0., prec));
	CHECK_THAT(virt_photon_from_virt_photon.t.y, AbsMatcher<Real>(0., prec));
	CHECK_THAT(virt_photon_from_virt_photon.t.z, AbsMatcher<Real>(0., prec));
	CHECK_THAT(virt_photon_from_virt_photon.x.t, AbsMatcher<Real>(0., prec));
	CHECK_THAT(virt_photon_from_virt_photon.x.x, AbsMatcher<Real>(-1., prec));
	CHECK_THAT(virt_photon_from_virt_photon.x.y, AbsMatcher<Real>(0., prec));
	CHECK_THAT(virt_photon_from_virt_photon.x.z, AbsMatcher<Real>(0., prec));
	CHECK_THAT(virt_photon_from_virt_photon.y.t, AbsMatcher<Real>(0., prec));
	CHECK_THAT(virt_photon_from_virt_photon.y.x, AbsMatcher<Real>(0., prec));
	CHECK_THAT(virt_photon_from_virt_photon.y.y, AbsMatcher<Real>(-1., prec));
	CHECK_THAT(virt_photon_from_virt_photon.y.z, AbsMatcher<Real>(0., prec));
	CHECK_THAT(virt_photon_from_virt_photon.z.t, AbsMatcher<Real>(0., prec));
	CHECK_THAT(virt_photon_from_virt_photon.z.x, AbsMatcher<Real>(0., prec));
	CHECK_THAT(virt_photon_from_virt_photon.z.y, AbsMatcher<Real>(0., prec));
	CHECK_THAT(virt_photon_from_virt_photon.z.z, AbsMatcher<Real>(-1., prec));

	// Check determinants.
	CHECK_THAT(target_from_lepton.det(), RelMatcher<Real>(1., prec));
	CHECK_THAT(target_from_hadron.det(), RelMatcher<Real>(1., prec));
	CHECK_THAT(target_from_virt_photon.det(), RelMatcher<Real>(1., prec));

	// Check individual frames.
	math::Vec4 q_test, k2_test, ph_test;
	q_test = target_from_lepton.transpose() * q;
	CHECK_THAT(q_test.t, RelMatcher<Real>(q.t, prec));
	CHECK_THAT(q_test.x, AbsMatcher<Real>(0., prec));
	CHECK_THAT(q_test.y, AbsMatcher<Real>(0., prec));
	CHECK_THAT(q_test.z, RelMatcher<Real>(q.r().norm(), prec));
	k2_test = target_from_lepton.transpose() * k2;
	CHECK_THAT(k2_test.t, RelMatcher<Real>(k2.t, prec));
	CHECK_THAT(k2_test.y, AbsMatcher<Real>(0., prec));

	q_test = target_from_hadron.transpose() * q;
	CHECK_THAT(q_test.t, RelMatcher<Real>(q.t, prec));
	CHECK_THAT(q_test.x, AbsMatcher<Real>(0., prec));
	CHECK_THAT(q_test.y, AbsMatcher<Real>(0., prec));
	CHECK_THAT(q_test.z, RelMatcher<Real>(q.r().norm(), prec));
	ph_test = target_from_hadron.transpose() * ph;
	CHECK_THAT(ph_test.t, RelMatcher<Real>(ph.t, prec));
	CHECK_THAT(ph_test.y, AbsMatcher<Real>(0., prec));

	q_test = target_from_virt_photon.transpose() * q;
	CHECK_THAT(q_test.t, AbsMatcher<Real>(0., prec));
	CHECK_THAT(q_test.x, AbsMatcher<Real>(0., prec));
	CHECK_THAT(q_test.y, AbsMatcher<Real>(0., prec));
	CHECK_THAT(q_test.z, RelMatcher<Real>(kin.Q, prec));
	ph_test = target_from_virt_photon.transpose() * ph;
	CHECK_THAT(ph_test.y, AbsMatcher<Real>(0., prec));
}

TEST_CASE(
		"Radiative reference frame checks",
		"[frame]") {
	InputRad input = GENERATE(
		from_stream<InputRad>(
			std::move(std::ifstream("data/phase_space_rad_vals.dat")),
			true));

	Real E_b = input.beam_energy;
	Real M = constant::MASS_P;
	constant::Lepton lep;
	if (input.particle_id == 'e') {
		lep = constant::Lepton::E;
	} else if (input.particle_id == 'm') {
		lep = constant::Lepton::MU;
	} else if (input.particle_id == 't') {
		lep = constant::Lepton::TAU;
	}
	Real M_th = constant::MASS_P + constant::MASS_PI_0;
	kin::Initial initial_state(constant::Nucleus::P, lep, E_b);
	kin::PhaseSpaceRad phase_space = input.phase_space;
	kin::KinematicsRad kin(initial_state, phase_space, constant::Hadron::PI_P, M_th);
	kin::FinalRad final_state(initial_state, math::Vec3::Y, kin);
	// Reference frames.
	math::Transform4 target_from_shift = frame::target_from_shift(kin);
	math::Transform4 target_from_hadron = frame::target_from_hadron(kin.project_shift());
	math::Transform4 target_from_real_photon = frame::target_from_real_photon(kin);
	// Get 4-momenta of particles.
	math::Vec4 p = initial_state.p;
	math::Vec4 k1 = initial_state.k1;
	math::Vec4 q = final_state.q;
	math::Vec4 k2 = final_state.k2;
	math::Vec4 k = final_state.k;
	math::Vec4 ph = final_state.ph;

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

	Real prec = 1e4 * std::numeric_limits<Real>::epsilon();

	// Check orthogonality.
	math::Transform4 shift_from_shift =
		target_from_shift.transpose() * target_from_shift;
	math::Transform4 real_photon_from_real_photon =
		target_from_real_photon.transpose() * target_from_real_photon;

	CHECK_THAT(shift_from_shift.t.t, AbsMatcher<Real>(1., prec));
	CHECK_THAT(shift_from_shift.t.x, AbsMatcher<Real>(0., prec));
	CHECK_THAT(shift_from_shift.t.y, AbsMatcher<Real>(0., prec));
	CHECK_THAT(shift_from_shift.t.z, AbsMatcher<Real>(0., prec));
	CHECK_THAT(shift_from_shift.x.t, AbsMatcher<Real>(0., prec));
	CHECK_THAT(shift_from_shift.x.x, AbsMatcher<Real>(-1., prec));
	CHECK_THAT(shift_from_shift.x.y, AbsMatcher<Real>(0., prec));
	CHECK_THAT(shift_from_shift.x.z, AbsMatcher<Real>(0., prec));
	CHECK_THAT(shift_from_shift.y.t, AbsMatcher<Real>(0., prec));
	CHECK_THAT(shift_from_shift.y.x, AbsMatcher<Real>(0., prec));
	CHECK_THAT(shift_from_shift.y.y, AbsMatcher<Real>(-1., prec));
	CHECK_THAT(shift_from_shift.y.z, AbsMatcher<Real>(0., prec));
	CHECK_THAT(shift_from_shift.z.t, AbsMatcher<Real>(0., prec));
	CHECK_THAT(shift_from_shift.z.x, AbsMatcher<Real>(0., prec));
	CHECK_THAT(shift_from_shift.z.y, AbsMatcher<Real>(0., prec));
	CHECK_THAT(shift_from_shift.z.z, AbsMatcher<Real>(-1., prec));

	CHECK_THAT(real_photon_from_real_photon.t.t, AbsMatcher<Real>(1., prec));
	CHECK_THAT(real_photon_from_real_photon.t.x, AbsMatcher<Real>(0., prec));
	CHECK_THAT(real_photon_from_real_photon.t.y, AbsMatcher<Real>(0., prec));
	CHECK_THAT(real_photon_from_real_photon.t.z, AbsMatcher<Real>(0., prec));
	CHECK_THAT(real_photon_from_real_photon.x.t, AbsMatcher<Real>(0., prec));
	CHECK_THAT(real_photon_from_real_photon.x.x, AbsMatcher<Real>(-1., prec));
	CHECK_THAT(real_photon_from_real_photon.x.y, AbsMatcher<Real>(0., prec));
	CHECK_THAT(real_photon_from_real_photon.x.z, AbsMatcher<Real>(0., prec));
	CHECK_THAT(real_photon_from_real_photon.y.t, AbsMatcher<Real>(0., prec));
	CHECK_THAT(real_photon_from_real_photon.y.x, AbsMatcher<Real>(0., prec));
	CHECK_THAT(real_photon_from_real_photon.y.y, AbsMatcher<Real>(-1., prec));
	CHECK_THAT(real_photon_from_real_photon.y.z, AbsMatcher<Real>(0., prec));
	CHECK_THAT(real_photon_from_real_photon.z.t, AbsMatcher<Real>(0., prec));
	CHECK_THAT(real_photon_from_real_photon.z.x, AbsMatcher<Real>(0., prec));
	CHECK_THAT(real_photon_from_real_photon.z.y, AbsMatcher<Real>(0., prec));
	CHECK_THAT(real_photon_from_real_photon.z.z, AbsMatcher<Real>(-1., prec));

	// Check determinants.
	CHECK_THAT(target_from_shift.det(), RelMatcher<Real>(1., prec));
	CHECK_THAT(target_from_real_photon.det(), RelMatcher<Real>(1., prec));

	// Check individual frames.
	math::Vec4 q_test, k_test;
	q_test = target_from_real_photon.transpose() * q;
	CHECK_THAT(q_test.t, RelMatcher<Real>(q.t, prec));
	CHECK_THAT(q_test.x, AbsMatcher<Real>(0., prec));
	CHECK_THAT(q_test.y, AbsMatcher<Real>(0., prec));
	CHECK_THAT(q_test.z, RelMatcher<Real>(q.r().norm(), prec));
	k_test = target_from_real_photon.transpose() * k;
	CHECK_THAT(k_test.t, RelMatcher<Real>(k.t, prec));
	CHECK_THAT(k_test.y, AbsMatcher<Real>(0., prec));

	// Check that the shift frame is the same as the one obtained by using
	// `project_shift` with the hadron frame.
	// TODO: Fix up this test so it can actually be used. Right now, it fails
	// because the precision requirements are too stringent since after
	// projection, `target_from_hadron` is no longer very accurate.
	/*
	CHECK_THAT(target_from_shift.t.t, RelMatcher<Real>(target_from_hadron.t.t, prec));
	CHECK_THAT(target_from_shift.t.x, RelMatcher<Real>(target_from_hadron.t.x, prec));
	CHECK_THAT(target_from_shift.t.y, RelMatcher<Real>(target_from_hadron.t.y, prec));
	CHECK_THAT(target_from_shift.t.z, RelMatcher<Real>(target_from_hadron.t.z, prec));
	CHECK_THAT(target_from_shift.x.t, RelMatcher<Real>(target_from_hadron.x.t, prec));
	CHECK_THAT(target_from_shift.x.x, RelMatcher<Real>(target_from_hadron.x.x, prec));
	CHECK_THAT(target_from_shift.x.y, RelMatcher<Real>(target_from_hadron.x.y, prec));
	CHECK_THAT(target_from_shift.x.z, RelMatcher<Real>(target_from_hadron.x.z, prec));
	CHECK_THAT(target_from_shift.y.t, RelMatcher<Real>(target_from_hadron.y.t, prec));
	CHECK_THAT(target_from_shift.y.x, RelMatcher<Real>(target_from_hadron.y.x, prec));
	CHECK_THAT(target_from_shift.y.y, RelMatcher<Real>(target_from_hadron.y.y, prec));
	CHECK_THAT(target_from_shift.y.z, RelMatcher<Real>(target_from_hadron.y.z, prec));
	CHECK_THAT(target_from_shift.z.t, RelMatcher<Real>(target_from_hadron.z.t, prec));
	CHECK_THAT(target_from_shift.z.x, RelMatcher<Real>(target_from_hadron.z.x, prec));
	CHECK_THAT(target_from_shift.z.y, RelMatcher<Real>(target_from_hadron.z.y, prec));
	CHECK_THAT(target_from_shift.z.z, RelMatcher<Real>(target_from_hadron.z.z, prec));
	*/
}

