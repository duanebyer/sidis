#include <catch2/catch.hpp>

#include <fstream>
#include <istream>
#include <limits>
#include <sstream>
#include <utility>

#include <sidis/constant.hpp>
#include <sidis/frame.hpp>
#include <sidis/kinematics.hpp>
#include <sidis/particle.hpp>
#include <sidis/transform.hpp>
#include <sidis/vector.hpp>
#include <sidis/extra/math.hpp>

#include "abs_matcher.hpp"
#include "phase_space_generator.hpp"
#include "rel_matcher.hpp"
#include "stream_generator.hpp"

using namespace sidis;

namespace {

struct Input {
	char particle_id;
	Real beam_energy;
	kin::PhaseSpace ph_space;
};

struct InputRad {
	char particle_id;
	Real beam_energy;
	kin::PhaseSpaceRad ph_space;
};

std::istream& operator>>(std::istream& in, Input& input) {
	in >> input.particle_id;
	if (!(input.particle_id == 'e'
			|| input.particle_id == 'm'
			|| input.particle_id == 't')) {
		in.setstate(std::ios_base::failbit);
	}
	in >> input.beam_energy;
	in >> input.ph_space.x;
	in >> input.ph_space.y;
	in >> input.ph_space.z;
	in >> input.ph_space.ph_t_sq;
	in >> input.ph_space.phi_h;
	in >> input.ph_space.phi;
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
	in >> input.ph_space.x;
	in >> input.ph_space.y;
	in >> input.ph_space.z;
	in >> input.ph_space.ph_t_sq;
	in >> input.ph_space.phi_h;
	in >> input.ph_space.phi;
	in >> input.ph_space.tau;
	in >> input.ph_space.R;
	in >> input.ph_space.phi_k;
	return in;
}

Real norm_euc(math::Vec4 vec) {
	return std::hypot(vec.t, vec.r().norm());
}

void test_kin_nrad(
		kin::Initial init,
		kin::Kinematics kin,
		bool complete,
		Real rel_prec=1e4) {
	kin::Final fin(init, math::VEC3_Y, kin);
	// Get 4-momenta of particles.
	math::Vec4 p = init.p;
	math::Vec4 k1 = init.k1;
	math::Vec4 q = fin.q;
	math::Vec4 k2 = fin.k2;
	math::Vec4 ph = fin.ph;
	math::Vec4 px = (p + k1) - (k2 + ph);
	// Basis vectors for angle checks.
	math::Vec3 e_y = cross(q.r(), k1.r()).unit();
	math::Vec3 e_x = cross(e_y, q.r()).unit();

	// Do comparisons.
	Real prec = rel_prec*std::numeric_limits<Real>::epsilon();

	// Kinematic variables.
	CHECK_THAT(
		-dot(q, q)/(2.*dot(q, p)),
		RelMatcher<Real>(
			kin.x,
			2.*prec*std::hypot(2.*q.t*q.t/dot(q, q), q.t*p.t/dot(q, p))));
	CHECK_THAT(
		dot(q, p)/dot(k1, p),
		RelMatcher<Real>(
			kin.y,
			2.*prec*std::hypot(q.t*p.t/dot(q, p), k1.t*p.t/dot(k1, p))));
	CHECK_THAT(
		dot(ph, p)/dot(p, q),
		RelMatcher<Real>(
			kin.z,
			2.*prec*std::hypot(ph.t*p.t/dot(ph, p), p.t*q.t/dot(p, q))));
	CHECK_THAT(
		(q - ph).norm_sq(),
		AbsMatcher<Real>(
			kin.t,
			6.*prec*std::hypot(q.t, ph.t)*norm_euc(q - ph)));
	CHECK_THAT(
		(q + p).norm_sq(),
		AbsMatcher<Real>(
			kin.W_sq,
			6.*prec*std::hypot(q.t, p.t)*norm_euc(q + p)));

	CHECK_THAT(
		2.*dot(p, k1),
		AbsMatcher<Real>(kin.S, 4.*prec*p.t*k1.t));
	CHECK_THAT(
		-q.norm_sq(),
		AbsMatcher<Real>(kin.Q_sq, 4.*prec*q.t*q.t));
	CHECK_THAT(
		2.*dot(p, k2),
		AbsMatcher<Real>(kin.X, 4.*prec*p.t*k2.t));
	CHECK_THAT(
		2.*dot(p, q),
		AbsMatcher<Real>(kin.S_x, 4.*prec*p.t*q.t));
	CHECK_THAT(
		2.*dot(k1, ph),
		AbsMatcher<Real>(kin.V_1, 4.*prec*k1.t*ph.t));
	CHECK_THAT(
		2.*dot(k2, ph),
		AbsMatcher<Real>(kin.V_2, 4.*prec*k2.t*ph.t));
	CHECK_THAT(
		dot(q, ph),
		AbsMatcher<Real>(kin.V_m, 2.*prec*q.t*ph.t));

	// 3-momenta magnitudes.
	CHECK_THAT(
		k1.r().norm(),
		RelMatcher<Real>(kin.lambda_S_sqrt/(2.*kin.M), 2.*prec));
	CHECK_THAT(
		q.r().norm(),
		RelMatcher<Real>(kin.lambda_Y_sqrt/(2.*kin.M), 2.*prec));

	// 3-momenta components.
	CHECK_THAT(
		dot(ph.r(), q.r().unit()),
		RelMatcher<Real>(kin.ph_l, 2.*prec));
	CHECK_THAT(
		cross(ph.r(), q.r().unit()).norm(),
		RelMatcher<Real>(kin.ph_t, 2.*prec));
	CHECK_THAT(
		cross(k1.r(), q.r().unit()).norm(),
		RelMatcher<Real>(kin.k1_t, 2.*prec));

	if (complete) {
		CHECK_THAT(
			cross(k2.r(), q.r().unit()).norm(),
			RelMatcher<Real>(kin.k1_t, 2.*prec));
	}
	CHECK_THAT(
		dot(q.r(), k1.r().unit()),
		RelMatcher<Real>(kin.q_l, 2.*prec));
	CHECK_THAT(
		cross(q.r(), k1.r().unit()).norm(),
		RelMatcher<Real>(kin.q_t, 2.*prec));

	// Volume parts.
	CHECK_THAT(
		dot(cross(k1.r(), q.r()), ph.r()),
		RelMatcher<Real>(kin.vol_phi_h/kin.M, 4.*prec));

	// Angles.
	CHECK_THAT(
		std::atan2(dot(e_y, ph.r()), dot(e_x, ph.r())),
		RelMatcher<Real>(kin.phi_h, 2.*prec));
	CHECK_THAT(
		std::atan2(k2.x, k2.y),
		RelMatcher<Real>(kin.phi, 2.*prec));
	CHECK_THAT(
		std::atan2(-q.x, q.y),
		RelMatcher<Real>(kin.phi_q, 2.*2.*prec));

	// Completeness.
	if (complete) {
		math::Vec4 q_test = k1 - k2;
		CHECK_THAT(
			q_test.t,
			AbsMatcher<Real>(q.t, prec*std::hypot(k1.t, k2.t)));
		CHECK_THAT(
			q_test.x,
			AbsMatcher<Real>(q.x, prec*std::hypot(k1.x, k2.x)));
		CHECK_THAT(
			q_test.y,
			AbsMatcher<Real>(q.y, prec*std::hypot(k1.y, k2.y)));
		CHECK_THAT(
			q_test.z,
			AbsMatcher<Real>(q.z, prec*std::hypot(k1.z, k2.z)));
	}

	// Conservation.
	math::Vec4 p_tot_i = p + k1;
	math::Vec4 p_tot_f = k2 + ph + px;
	CHECK_THAT(
		p_tot_i.t,
		AbsMatcher<Real>(p_tot_f.t, prec*std::sqrt(p.t*p.t + k1.t*k1.t + k2.t*k2.t + ph.t*ph.t + px.t*px.t)));
	CHECK_THAT(
		p_tot_i.x,
		AbsMatcher<Real>(p_tot_f.x, prec*std::sqrt(p.x*p.x + k1.x*k1.x + k2.x*k2.x + ph.x*ph.x + px.x*px.x)));
	CHECK_THAT(
		p_tot_i.y,
		AbsMatcher<Real>(p_tot_f.y, prec*std::sqrt(p.y*p.y + k1.y*k1.y + k2.y*k2.y + ph.y*ph.y + px.y*px.y)));
	CHECK_THAT(
		p_tot_i.z,
		AbsMatcher<Real>(p_tot_f.z, prec*std::sqrt(p.z*p.z + k1.z*k1.z + k2.z*k2.z + ph.z*ph.z + px.z*px.z)));

	// Particle masses.
	CHECK_THAT(
		p.norm(),
		AbsMatcher<Real>(kin.M, 4.*prec*p.t*p.t/kin.M));
	CHECK_THAT(
		k1.norm(),
		AbsMatcher<Real>(kin.m, 4.*prec*k1.t*k1.t/kin.m));
	CHECK_THAT(
		k2.norm(),
		AbsMatcher<Real>(kin.m, 4.*prec*k2.t*k2.t/kin.m));
	CHECK_THAT(
		ph.norm(),
		AbsMatcher<Real>(kin.mh, 4.*prec*ph.t*ph.t/kin.mh));
	if (complete) {
		CHECK_THAT(
			px.norm(),
			AbsMatcher<Real>(kin.mx, 4.*prec*px.t*px.t/kin.mx));
	}
}

void test_kin_rad(
		kin::Initial init,
		kin::KinematicsRad kin,
		bool complete,
		Real rel_prec=1e4) {
	kin::FinalRad fin(init, math::VEC3_ZERO, kin);
	// Get 4-momenta of particles.
	math::Vec4 p = init.p;
	math::Vec4 k1 = init.k1;
	math::Vec4 q = fin.q;
	math::Vec4 k2 = fin.k2;
	math::Vec4 k = fin.k;
	math::Vec4 ph = fin.ph;
	math::Vec4 px = (p + k1) - (k2 + ph);
	// Basis vectors for angle checks.
	math::Vec3 e_y = cross(k1.r(), k2.r()).unit();
	math::Vec3 e_x = cross(e_y, q.r()).unit();
	math::Vec3 shift_e_y = cross(k1.r(), (k2 + k).r()).unit();
	math::Vec3 shift_e_x = cross(shift_e_y, (q - k).r()).unit();

	// Do comparisons.
	Real prec = rel_prec*std::numeric_limits<Real>::epsilon();

	// Kinematic variables.
	CHECK_THAT(
		dot(k, q)/dot(k, p),
		RelMatcher<Real>(
			kin.tau,
			2.*prec*std::hypot(k.t*q.t/dot(k, q), k.t*p.t/dot(k, p))));
	CHECK_THAT(
		2.*dot(k, p),
		AbsMatcher<Real>(kin.R, 4.*prec*k.t*p.t));
	CHECK_THAT(
		dot(k1, k)/dot(p, k),
		RelMatcher<Real>(
			kin.z_1,
			2.*prec*std::hypot(k1.t*k.t/dot(k1, k), p.t*k.t/dot(p, k))));
	CHECK_THAT(
		dot(k2, k)/dot(p, k),
		RelMatcher<Real>(
			kin.z_2,
			2.*prec*std::hypot(k2.t*k.t/dot(k2, k), p.t*k.t/dot(p, k))));
	CHECK_THAT(
		dot(k, ph)/dot(k, p),
		RelMatcher<Real>(
			kin.mu,
			2.*prec*std::hypot(k.t*ph.t/dot(k, ph), k.t*p.t/dot(k, p))));
	CHECK_THAT(
		(k1 / kin.z_1 - k2 / kin.z_2).norm_sq(),
		AbsMatcher<Real>(
			kin.F_IR,
			6.*prec*std::hypot(k1.t / kin.z_1, k2.t / kin.z_2)*norm_euc(k1 / kin.z_1 - k2 / kin.z_2)));

	// 3-momenta dot products.
	CHECK_THAT(
		dot(k.r(), ph.r()),
		RelMatcher<Real>(kin.lambda_RV/(4.*kin.M*kin.M), 2.*prec));
	CHECK_THAT(
		dot(k.r(), q.r()),
		RelMatcher<Real>(kin.lambda_RY/(4.*kin.M*kin.M), 2.*prec));
	CHECK_THAT(
		dot(q.r(), ph.r()),
		RelMatcher<Real>(kin.lambda_V/(4.*kin.M*kin.M), 2.*prec));
	CHECK_THAT(
		ph.r().norm_sq(),
		RelMatcher<Real>(kin.lambda_H/(4.*kin.M*kin.M), 2.*prec));

	// Volume parts.
	CHECK_THAT(
		dot(cross(k1.r(), q.r()), k.r()),
		RelMatcher<Real>(kin.R*kin.vol_phi_k_R/kin.M, 4.*prec));
	CHECK_THAT(
		dot(cross(k.r(), ph.r()), q.r()),
		RelMatcher<Real>(kin.vol_phi_hk/kin.M, 4.*prec));

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
		RelMatcher<Real>(
			kin.shift_x,
			2.*prec*std::hypot(2.*math::sq(q.t - k.t)/dot(q - k, q - k), (q.t - k.t)*p.t/dot(q - k, p))));
	CHECK_THAT(
		dot(ph, p)/dot(p, q - k),
		RelMatcher<Real>(
			kin.shift_z,
			2.*prec*std::hypot(ph.t*p.t/dot(ph, p), p.t*(q.t - k.t)/dot(p, q - k))));
	CHECK_THAT(
		(q - k - ph).norm_sq(),
		RelMatcher<Real>(
			kin.shift_t,
			6.*prec*std::hypot(q.t - k.t, ph.t)*norm_euc(q - k - ph)));
	CHECK_THAT(
		-(q - k).norm_sq(),
		AbsMatcher<Real>(kin.shift_Q_sq, 4.*prec*math::sq(q.t - k.t)));
	CHECK_THAT(
		2.*dot(p, q - k),
		AbsMatcher<Real>(kin.shift_S_x, 4.*prec*p.t*(q.t - k.t)));
	CHECK_THAT(
		dot(q - k, ph),
		AbsMatcher<Real>(kin.shift_V_m, 2.*prec*(q.t - k.t)*ph.t));

	// Shifted 3-momenta magnitudes.
	CHECK_THAT(
		(q - k).r().norm(),
		RelMatcher<Real>(kin.shift_lambda_Y_sqrt/(2.*kin.M), 2.*prec));

	// Shifted 3-momenta components.
	CHECK_THAT(
		dot(ph.r(), (q - k).r().unit()),
		RelMatcher<Real>(kin.shift_ph_l, 2.*prec));
	CHECK_THAT(
		cross(ph.r(), (q - k).r().unit()).norm(),
		RelMatcher<Real>(kin.shift_ph_t, 2.*prec));
	CHECK_THAT(
		cross(k1.r(), (q - k).r().unit()).norm(),
		RelMatcher<Real>(kin.shift_k1_t, 2.*prec));
	CHECK_THAT(
		dot((q - k).r(), k1.r().unit()),
		RelMatcher<Real>(kin.shift_q_l, 2.*prec));
	CHECK_THAT(
		cross((q - k).r(), k1.r().unit()).norm(),
		RelMatcher<Real>(kin.shift_q_t, 2.*prec));

	// Shifted volume parts.
	CHECK_THAT(
		dot(cross(k1.r(), (q - k).r()), ph.r()),
		RelMatcher<Real>(kin.shift_vol_phi_h/kin.M, 4.*prec));

	// Shifted angles.
	CHECK_THAT(
		std::atan2(dot(shift_e_y, ph.r()), dot(shift_e_x, ph.r())),
		RelMatcher<Real>(kin.shift_phi_h, 2.*prec));
	CHECK_THAT(
		std::atan2((q - k).x, -(q - k).y),
		RelMatcher<Real>(kin.shift_phi_q, 2.*2.*prec));
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
	Real M = MASS_P;
	part::Lepton lep;
	if (input.particle_id == 'e') {
		lep = part::Lepton::E;
	} else if (input.particle_id == 'm') {
		lep = part::Lepton::MU;
	} else if (input.particle_id == 't') {
		lep = part::Lepton::TAU;
	}
	Real Mth = MASS_P + MASS_PI_0;
	part::Particles ps(part::Nucleus::P, lep, part::Hadron::PI_P, Mth);
	Real S = 2.*ps.M*E_b;
	kin::Initial init(ps, E_b);
	kin::PhaseSpace ph_space = input.ph_space;
	kin::Kinematics kin(ps, S, ph_space);

	// Print state information.
	std::stringstream ss;
	ss
		<< "pid   = " << input.particle_id   << std::endl
		<< "E_b   = " << E_b                 << std::endl
		<< "x     = " << ph_space.x       << std::endl
		<< "y     = " << ph_space.y       << std::endl
		<< "z     = " << ph_space.z       << std::endl
		<< "ph_t² = " << ph_space.ph_t_sq << std::endl
		<< "φ_h   = " << ph_space.phi_h   << std::endl
		<< "φ     = " << ph_space.phi     << std::endl;
	INFO(ss.str());

	test_kin_nrad(init, kin, true);
}

TEST_CASE(
		"Non-radiative shifted kinematics checks",
		"[kin]") {
	InputRad input = GENERATE(
		from_stream<InputRad>(
			std::move(std::ifstream("data/phase_space_rad_vals.dat")),
			true));
	Real E_b = input.beam_energy;
	Real M = MASS_P;
	part::Lepton lep;
	if (input.particle_id == 'e') {
		lep = part::Lepton::E;
	} else if (input.particle_id == 'm') {
		lep = part::Lepton::MU;
	} else if (input.particle_id == 't') {
		lep = part::Lepton::TAU;
	}
	Real Mth = MASS_P + MASS_PI_0;
	part::Particles ps(part::Nucleus::P, lep, part::Hadron::PI_P, Mth);
	Real S = 2.*ps.M*E_b;
	kin::Initial init(ps, E_b);
	kin::PhaseSpaceRad ph_space = input.ph_space;
	kin::KinematicsRad kin(ps, S, ph_space);

	// Print state information.
	std::stringstream ss;
	ss
		<< "pid   = " << input.particle_id   << std::endl
		<< "E_b   = " << E_b                 << std::endl
		<< "x     = " << ph_space.x       << std::endl
		<< "y     = " << ph_space.y       << std::endl
		<< "z     = " << ph_space.z       << std::endl
		<< "ph_t² = " << ph_space.ph_t_sq << std::endl
		<< "φ_h   = " << ph_space.phi_h   << std::endl
		<< "φ     = " << ph_space.phi     << std::endl
		<< "τ     = " << ph_space.tau     << std::endl
		<< "R     = " << ph_space.R       << std::endl
		<< "φ_k   = " << ph_space.phi_k   << std::endl;
	INFO(ss.str());

	test_kin_nrad(init, kin.project_shift(), false);
}

TEST_CASE(
		"Radiative kinematics checks",
		"[kin]") {
	InputRad input = GENERATE(
		from_stream<InputRad>(
			std::move(std::ifstream("data/phase_space_rad_vals.dat")),
			true));
	Real E_b = input.beam_energy;
	Real M = MASS_P;
	part::Lepton lep;
	if (input.particle_id == 'e') {
		lep = part::Lepton::E;
	} else if (input.particle_id == 'm') {
		lep = part::Lepton::MU;
	} else if (input.particle_id == 't') {
		lep = part::Lepton::TAU;
	}
	Real Mth = MASS_P + MASS_PI_0;
	part::Particles ps(part::Nucleus::P, lep, part::Hadron::PI_P, Mth);
	Real S = 2.*ps.M*E_b;
	kin::Initial init(ps, E_b);
	kin::PhaseSpaceRad ph_space = input.ph_space;
	kin::KinematicsRad kin(ps, S, ph_space);

	// Print state information.
	std::stringstream ss;
	ss
		<< "pid   = " << input.particle_id   << std::endl
		<< "E_b   = " << E_b                 << std::endl
		<< "x     = " << ph_space.x       << std::endl
		<< "y     = " << ph_space.y       << std::endl
		<< "z     = " << ph_space.z       << std::endl
		<< "ph_t² = " << ph_space.ph_t_sq << std::endl
		<< "φ_h   = " << ph_space.phi_h   << std::endl
		<< "φ     = " << ph_space.phi     << std::endl
		<< "τ     = " << ph_space.tau     << std::endl
		<< "R     = " << ph_space.R       << std::endl
		<< "φ_k   = " << ph_space.phi_k   << std::endl;
	INFO(ss.str());

	test_kin_rad(init, kin, true, 1e4);
}

// The following `[kin-rand]` tests will likely have a handful of failures,
// corresponding to tested points that have bad numerical properties (resulting
// in larger than expected inaccuracies in the result). Because of this, they
// are disabled by default, but they can still be useful in verifying that
// nothing really bad is going with kinematic calculations.
TEST_CASE(
		"Random phase space bounds outer check",
		"[.][kin-rand]") {
	// Generate a random point just outside the phase space boundary, then check
	// that it is not kinematically valid.
	Real E_b = GENERATE(3., 12., 140., 3000., 12000.);
	Real Mth = MASS_P + MASS_PI_0;
	part::Nucleus target = part::Nucleus::P;
	part::Lepton lepton = part::Lepton::TAU;
	part::Hadron hadron = part::Hadron::PI_P;
	part::Particles ps(target, lepton, hadron, Mth);
	Real S = 2.*ps.M*E_b;
	kin::Initial init(ps, E_b);
	kin::PhaseSpace ph_space = GENERATE_COPY(
		take(200000, gen_phase_space_surface(ps, S, -0.0001)));
	kin::Kinematics kin(ps, S, ph_space);

	std::stringstream ss;
	ss
		<< "E_b   = " << E_b                 << std::endl
		<< "x     = " << ph_space.x       << std::endl
		<< "y     = " << ph_space.y       << std::endl
		<< "z     = " << ph_space.z       << std::endl
		<< "ph_t² = " << ph_space.ph_t_sq << std::endl
		<< "φ_h   = " << ph_space.phi_h   << std::endl
		<< "φ     = " << ph_space.phi     << std::endl;
	INFO(ss.str());

	CHECK(!cut::valid(kin));
}

TEST_CASE(
		"Random phase space bounds inner check",
		"[.][kin-rand]") {
	// Generate a random point just inside phase space boundary, then check that
	// it is kinematically valid.
	Real E_b = GENERATE(3., 12., 140., 3000., 12000.);
	Real Mth = MASS_P + MASS_PI_0;
	part::Nucleus target = part::Nucleus::P;
	part::Lepton lepton = part::Lepton::TAU;
	part::Hadron hadron = part::Hadron::PI_P;
	part::Particles ps(target, lepton, hadron, Mth);
	Real S = 2.*ps.M*E_b;
	kin::Initial init(ps, E_b);
	kin::PhaseSpace ph_space = GENERATE_COPY(
		take(200000, gen_phase_space_surface(ps, S, 0.0001)));
	kin::Kinematics kin(ps, S, ph_space);

	std::stringstream ss;
	ss
		<< "E_b   = " << E_b                 << std::endl
		<< "x     = " << ph_space.x       << std::endl
		<< "y     = " << ph_space.y       << std::endl
		<< "z     = " << ph_space.z       << std::endl
		<< "ph_t² = " << ph_space.ph_t_sq << std::endl
		<< "φ_h   = " << ph_space.phi_h   << std::endl
		<< "φ     = " << ph_space.phi     << std::endl;
	INFO(ss.str());

	REQUIRE(cut::valid(kin));
	test_kin_nrad(init, kin, true, 1e6);
}

TEST_CASE(
		"Random phase space points check",
		"[.][kin-rand]") {
	// We choose these conditions to ensure that the phase space is explored
	// even in the region between non-relativistic and ultra-relativistic.
	Real E_b = GENERATE(3., 12., 140., 3000., 12000.);
	Real Mth = MASS_P + MASS_PI_0;
	part::Nucleus target = part::Nucleus::P;
	part::Lepton lepton = part::Lepton::TAU;
	part::Hadron hadron = part::Hadron::PI_P;
	part::Particles ps(target, lepton, hadron, Mth);
	Real S = 2.*ps.M*E_b;
	kin::Initial init(ps, E_b);
	kin::PhaseSpace ph_space = GENERATE_COPY(
		take(200000, gen_phase_space(ps, S)));
	kin::Kinematics kin(ps, S, ph_space);

	std::stringstream ss;
	ss
		<< "pid   = " << name(lepton)        << std::endl
		<< "E_b   = " << E_b                 << std::endl
		<< "x     = " << ph_space.x       << std::endl
		<< "y     = " << ph_space.y       << std::endl
		<< "z     = " << ph_space.z       << std::endl
		<< "ph_t² = " << ph_space.ph_t_sq << std::endl
		<< "φ_h   = " << ph_space.phi_h   << std::endl
		<< "φ     = " << ph_space.phi     << std::endl;
	INFO(ss.str());

	test_kin_nrad(init, kin, true, 1e6);
}

TEST_CASE(
		"Target to lab frame checks",
		"[frame]") {
	math::Vec3 p(1.2, -0.5, 2.3);
	math::Vec3 k1(0.4, -2.2, 0.2);
	math::Vec3 pol(0.2, 0.1, -0.3);
	part::Particles ps(
		part::Nucleus::P,
		part::Lepton::MU,
		part::Hadron::PI_P,
		MASS_P + MASS_PI_0);
	kin::Initial init(ps, p, k1);

	// Construct the frames.
	math::Transform4 target_from_lab = frame::target_from_lab(init, pol);
	math::Transform4 lab_from_target = frame::lab_from_target(init, pol);
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
	math::Vec4 p_target = target_from_lab * init.p;
	math::Vec4 k1_target = target_from_lab * init.k1;
	CHECK_THAT(p_target.t, RelMatcher<Real>(MASS_P, prec));
	CHECK_THAT(p_target.x, AbsMatcher<Real>(0., prec));
	CHECK_THAT(p_target.y, AbsMatcher<Real>(0., prec));
	CHECK_THAT(p_target.z, AbsMatcher<Real>(0., prec));
	CHECK(k1_target.t > 0.);
	CHECK_THAT(k1_target.x, AbsMatcher<Real>(0., prec));
	CHECK_THAT(k1_target.y, AbsMatcher<Real>(0., prec));
	CHECK(k1_target.z > 0.);
	CHECK_THAT(k1_target.norm(), RelMatcher<Real>(MASS_MU, prec));

	// Check transforming polarization correctly. Remember that polarization is
	// given in the proton rest frame, so first have to boost it into proton
	// frame.
	math::Transform4 pol_bv(
		0., 0., 0., 0.,
		0., 0., -pol.z, pol.y,
		0., pol.z, 0., -pol.x,
		0., -pol.y, pol.x, 0.);
	math::Transform4 boost = math::Transform4::boost_to(init.p);
	math::Vec4 k1_boost = boost.transpose() * init.k1;
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
	part::Particles ps(
		part::Nucleus::P,
		part::Lepton::MU,
		part::Hadron::PI_P,
		MASS_P + MASS_PI_0);
	kin::Initial init(ps, 8.2);

	// Construct the frames.
	math::Vec3 pol = GENERATE(
		math::Vec3(0., 0., 0.),
		math::Vec3(0., 0., 1e-8),
		math::Vec3(0., 1e-9, 1e-8),
		math::Vec3(0.1, -0.2, 0.3),
		math::Vec3(0., 0., -0.4));
	math::Transform4 target_from_lab = frame::target_from_lab(init, pol);
	math::Transform4 lab_from_target = frame::lab_from_target(init, pol);
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
	Real M = MASS_P;
	part::Lepton lep;
	if (input.particle_id == 'e') {
		lep = part::Lepton::E;
	} else if (input.particle_id == 'm') {
		lep = part::Lepton::MU;
	} else if (input.particle_id == 't') {
		lep = part::Lepton::TAU;
	}
	Real Mth = MASS_P + MASS_PI_0;
	part::Particles ps(part::Nucleus::P, lep, part::Hadron::PI_P, Mth);
	Real S = 2.*ps.M*E_b;
	kin::Initial init(ps, E_b);
	kin::PhaseSpace ph_space = input.ph_space;
	kin::Kinematics kin(ps, S, ph_space);
	kin::Final fin(init, math::VEC3_Y, kin);
	// Reference frames.
	math::Transform4 target_from_lepton = frame::target_from_lepton(kin);
	math::Transform4 target_from_hadron = frame::target_from_hadron(kin);
	math::Transform4 target_from_virt_photon = frame::target_from_virt_photon(kin);
	// Get 4-momenta of particles.
	math::Vec4 p = init.p;
	math::Vec4 k1 = init.k1;
	math::Vec4 q = fin.q;
	math::Vec4 k2 = fin.k2;
	math::Vec4 ph = fin.ph;

	// Print state information.
	std::stringstream ss;
	ss
		<< "pid   = " << input.particle_id   << std::endl
		<< "E_b   = " << E_b                 << std::endl
		<< "x     = " << ph_space.x       << std::endl
		<< "y     = " << ph_space.y       << std::endl
		<< "z     = " << ph_space.z       << std::endl
		<< "ph_t² = " << ph_space.ph_t_sq << std::endl
		<< "φ_h   = " << ph_space.phi_h   << std::endl
		<< "φ     = " << ph_space.phi     << std::endl;
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
	Real M = MASS_P;
	part::Lepton lep;
	if (input.particle_id == 'e') {
		lep = part::Lepton::E;
	} else if (input.particle_id == 'm') {
		lep = part::Lepton::MU;
	} else if (input.particle_id == 't') {
		lep = part::Lepton::TAU;
	}
	Real Mth = MASS_P + MASS_PI_0;
	part::Particles ps(part::Nucleus::P, lep, part::Hadron::PI_P, Mth);
	Real S = 2.*ps.M*E_b;
	kin::Initial init(ps, E_b);
	kin::PhaseSpaceRad ph_space = input.ph_space;
	kin::KinematicsRad kin(ps, S, ph_space);
	kin::FinalRad fin(init, math::VEC3_Y, kin);
	// Reference frames.
	math::Transform4 target_from_shift = frame::target_from_shift(kin);
	math::Transform4 target_from_hadron = frame::target_from_hadron(kin.project_shift());
	math::Transform4 target_from_real_photon = frame::target_from_real_photon(kin);
	// Get 4-momenta of particles.
	math::Vec4 p = init.p;
	math::Vec4 k1 = init.k1;
	math::Vec4 q = fin.q;
	math::Vec4 k2 = fin.k2;
	math::Vec4 k = fin.k;
	math::Vec4 ph = fin.ph;

	// Print state information.
	std::stringstream ss;
	ss
		<< "pid   = " << input.particle_id   << std::endl
		<< "E_b   = " << E_b                 << std::endl
		<< "x     = " << ph_space.x       << std::endl
		<< "y     = " << ph_space.y       << std::endl
		<< "z     = " << ph_space.z       << std::endl
		<< "ph_t² = " << ph_space.ph_t_sq << std::endl
		<< "φ_h   = " << ph_space.phi_h   << std::endl
		<< "φ     = " << ph_space.phi     << std::endl
		<< "τ     = " << ph_space.tau     << std::endl
		<< "R     = " << ph_space.R       << std::endl
		<< "φ_k   = " << ph_space.phi_k   << std::endl;
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
}

