#include "sidis/cut.hpp"

#include <algorithm>
#include <cmath>

#include "sidis/constant.hpp"
#include "sidis/kinematics.hpp"
#include "sidis/particle.hpp"
#include "sidis/extra/math.hpp"

using namespace sidis;
using namespace sidis::cut;
using namespace sidis::kin;
using namespace sidis::math;
using namespace sidis::part;

namespace {

// To improve quality of numerical integration of the cross-section over the
// phase space, we introduce some possible variable transformations. These
// transformations are designed to "flatten" the integrand. They get used in the
// `take` methods below.

// Some constants used for the transformations.

// This is from emperical measurements. It doesn't have to be exactly right to
// still improve the integration substantially. Within an order or so will do.
static Real const PH_T_SQ_DECAY_WIDTH = 0.25;
// This x cutoff is chosen specifically because a lot of parameterizations of
// structure functions don't extend much below 1e-3.
static Real const X_CUTOFF = 1e-3;

// Trivial linear transformation.
struct Linear {
	Real u(Real x) const {
		return x;
	}
	Real x(Real u, Real* jac) const {
		*jac = 1.;
		return u;
	}
};

// Transforms according to reciprocal function.
struct Inverse {
	Real offset;
	Inverse(Real offset=0.) : offset(offset) { }
	Real u(Real x) const {
		return -1. / (x + offset);
	}
	Real x(Real u, Real* jac) const {
		*jac = sq(1. / u);
		return -1. / u - offset;
	}
};

// Transforms according to logarithm.
struct Log {
	Real offset;
	Log(Real offset=0.) : offset(offset) { }
	Real u(Real x) const {
		return std::log(x + offset);
	}
	Real x(Real u, Real* jac) const {
		Real exp = std::exp(u);
		*jac = exp;
		return exp - offset;
	}
};

// Transforms according to decaying exponential.
struct Decay {
	Real length;
	Decay(Real length=1.) : length(length) { }
	Real u(Real x) const {
		return -std::expm1(-x / length);
	}
	Real x(Real u, Real* jac) const {
		*jac = length / (1. - u);
		return -length * std::log1p(-u);
	}
};

// Transforms according to sigmoid. Useful for peaked functions.
struct Sigmoid {
	Real center;
	Real width;
	Sigmoid(Real center=0., Real width=1.) : center(center), width(width) { }
	Real u(Real x) const {
		// Use the `arcsinh` function here, as it has logarithmic tails which
		// allows for the non-peak part of the integrand to still be captured
		// accurately, even with fat tails.
		return std::asinh((x - center) / width);
	}
	Real x(Real u, Real* jac) const {
		*jac = width * std::cosh(u);
		return width * std::sinh(u) + center;
	}
};

// Transforms according to a double-sigmoid. Useful for doubly-peaked
// integrands. One of the peaks is specified in the regular way, and then the
// other peak is specified by the ratio of its center and width to those of the
// first peak.
struct Sigmoid2 {
	Real center;
	Real width;
	// Ratio between second sigmoid center and first one.
	Real center_r;
	// Ratio between second sigmoid width and first one.
	Real width_r;
	Sigmoid2(Real center=0., Real width=1., Real center_r=-1., Real width_r=1.) :
		center(center),
		width(width),
		center_r(center_r),
		width_r(width_r) { }
	Real u(Real x) const {
		Real t1 = (x - center) / width;
		Real t2 = (x - center * center_r) / (width * width_r);
		return std::asinh(t1) + std::asinh(t2);
	}
	Real x(Real u, Real* jac) const {
		// Against all odds, this function really can be inverted.
		Real sinhu = std::sinh(u);
		Real coshu = std::cosh(u);
		Real a = center_r + sq(width_r) + width_r * (1. + center_r) * coshu;
		Real b = 1. + sq(width_r) + 2. * width_r * coshu;
		Real c = sq((center_r - 1.) * (center / width));
		Real det = std::sqrt(b + c);
		Real a_p = width_r * (1. + center_r) * sinhu;
		Real b_p = 2. * width_r * sinhu;
		Real num = center * a + width * width_r * det * sinhu;
		*jac = (center * a_p
			+ width * width_r * (det * coshu + 0.5 * sinhu * b_p / det)
			- num * b_p / b) / b;
		return num / b;
	}
};

// Uses the specified transform to take a point `p` in the range `[0, 1)` and
// map it into another `bound`, with a resulting Jacobian `jac`.
template<typename T>
Real transform(T const& t, Real p, Bound bound, Real* jac) {
	Bound u_bound(t.u(bound.min()), t.u(bound.max()));
	Real u = u_bound.lerp(p);
	Real x = t.x(u, jac);
	*jac *= u_bound.size();
	return x;
}

}

Cut::Cut() :
	x(Bound::INVALID),
	y(Bound::INVALID),
	z(Bound::INVALID),
	ph_t_sq(Bound::INVALID),
	phi_h(Bound::INVALID),
	phi(Bound::INVALID),
	Q_sq(Bound::INVALID),
	t(Bound::INVALID),
	W_sq(Bound::INVALID),
	r(Bound::INVALID),
	mx_sq(Bound::INVALID),
	q_0(Bound::INVALID),
	k2_0(Bound::INVALID),
	ph_0(Bound::INVALID),
	theta_q(Bound::INVALID),
	theta_k2(Bound::INVALID),
	theta_h(Bound::INVALID) { }

CutRad::CutRad() :
	tau(Bound::INVALID),
	phi_k(Bound::INVALID),
	R(Bound::INVALID),
	k_0_bar(Bound::INVALID),
	k_0(Bound::INVALID),
	theta_k(Bound::INVALID) { }

// Bound of kinematic variables.
// TODO: Some of the calculations in this section are redundant with earlier
// kinematic calculations. This should be refactored to avoid that later.
Real cut::S_min(Particles const& ps) {
	return sq(ps.Mth + ps.mh) + 2.*ps.m*(ps.Mth + ps.mh) - sq(ps.M);
}
Bound cut::x_bound(Particles const& ps, Real S) {
	Real M = ps.M;
	Real m = ps.m;
	Real mh = ps.mh;
	Real Mth = ps.Mth;
	Real lambda_S = sq(S) - 4.*sq(M)*sq(m);
	Real L = S - (sq(Mth + mh) - sq(M));

	// Combined bound `mx >= Mth` and `sq(q_t) >= 0`.
	Real denom = 2.*(lambda_S + sq(M)*(sq(Mth + mh) - sq(M)));
	Real a = lambda_S - S*(sq(Mth + mh) - sq(M));
	Real b = std::sqrt(lambda_S*(L - 2.*m*(Mth + mh))*(L + 2.*m*(Mth + mh)));
	Real x_1 = (a - b)/denom;
	Real x_2 = (a + b)/denom;
	Bound kin_b(x_1, x_2);

	return Bound::UNIT & kin_b;
}
Bound cut::y_bound(Particles const& ps, Real S, Real x) {
	Real M = ps.M;
	Real m = ps.m;
	Real mh = ps.mh;
	Real Mth = ps.Mth;
	Real lambda_S = sq(S) - 4.*sq(M)*sq(m);

	// Bound `mx >= Mth`.
	Real px_threshold = (sq(Mth + mh) - sq(M))/((1. - x)*S);
	// Bound `sq(q_t) >= 0`.
	Real q_threshold = (x*lambda_S)/(S*(sq(M*x) + S*x + sq(m)));
	Bound kin_b(px_threshold, q_threshold);

	return Bound::UNIT & kin_b;
}
Bound cut::z_bound(Particles const& ps, Real S, Real x, Real y) {
	Real M = ps.M;
	Real mh = ps.mh;
	Real Mth = ps.Mth;
	Real S_x = S*y;
	Real Q_sq = S*x*y;
	Real lambda_Y = sq(S_x) + 4.*sq(M)*Q_sq;
	Real L = (1. - x)*S_x + sq(M);

	// Bound `mx >= Mth`.
	Real denom = 2.*S_x*L;
	Real a = (S_x + 2.*sq(M))*(L - sq(Mth) + sq(mh));
	Real b = std::sqrt(lambda_Y*(L - sq(Mth - mh))*(L - sq(Mth + mh)));
	Real z_crossover = 2.*sq(M)*(L - sq(Mth) + sq(mh))/(S_x*(S_x + 2.*sq(M)));
	Real z_0 = (2.*M*mh)/S_x;
	Real z_1 = (a - b)/denom;
	Real z_2 = (a + b)/denom;
	Real px_threshold_min = z_crossover > z_0 ? z_0 : z_1;
	Real px_threshold_max = z_2;
	Bound kin_b(px_threshold_min, px_threshold_max);

	return Bound::UNIT & kin_b;
}
Bound cut::ph_t_sq_bound(Particles const& ps, Real S, Real x, Real y, Real z) {
	Real M = ps.M;
	Real mh = ps.mh;
	Real Mth = ps.Mth;
	Real S_x = S*y;
	Real Q_sq = S*x*y;
	Real lambda_Y = sq(S_x) + 4.*sq(M)*Q_sq;
	Real L = (1. - x)*S_x + sq(M);
	Real ph_0 = (z*S_x)/(2.*M);

	// Bound `mx >= Mth`.
	Real det = (S_x + 2.*sq(M))*S_x/(2.*sq(M))*z - L + sq(Mth) - sq(mh);
	if (det < 0.) {
		det = 0.;
	}
	Real px_threshold = sq(ph_0) - sq(mh) - sq(M)/lambda_Y*sq(det);
	Bound kin_b(0., px_threshold);

	return Bound::POSITIVE & kin_b;
}
Bound cut::tau_bound(Kinematics const& kin) {
	// Equation [1.44].
	Real min = (kin.S_x - kin.lambda_Y_sqrt)/(2.*sq(kin.M));
	Real max = (kin.S_x + kin.lambda_Y_sqrt)/(2.*sq(kin.M));
	Bound kin_b(min, max);

	return Bound::FULL & kin_b;
}
Bound cut::R_bound(Kinematics const& kin, Real tau, Real phi_k) {
	Bound tau_b = tau_bound(kin);
	// Copied from kinematic calculations.
	Real mu = kin.ph_0/kin.M
		+ 1./kin.lambda_Y_sqrt*(
			(2.*tau*sq(kin.M) - kin.S_x)*kin.ph_l/kin.M
			- 2.*kin.M*kin.ph_t*std::cos(kin.phi_h - phi_k)*std::sqrt(
				(tau - tau_b.min())*(tau_b.max() - tau)));
	// Equation [1.44].
	Bound rad(0., 1./(1. + tau - mu)*(kin.mx_sq - sq(kin.Mth)));
	// It's unclear whether these bounds are necessary. For now, they are
	// commented out.
	/*
	// These additional bounds ensure that the shifted structure functions are
	// evaluated at valid points.
	// Shifted x and z variables must be between zero and one.
	Bound x(0., (1. - kin.x) / (1. + tau) * kin.S_x);
	Bound z(0., (1. - kin.z) * kin.S_x);
	// Shifted Q-squared must be positive.
	Bound Q_sq = tau < 0. ? Bound(0., -kin.Q_sq / tau) : Bound::POSITIVE;
	// TODO: Missing an additional bound to ensure shifted `ph_t_sq` is
	// positive.
	*/
	Bound kin_b = rad;
	return kin_b;
}

// TODO: Verify the extra kinematics cuts.
Bound cut::x_bound(Cut const& cut, Particles const& ps, Real S) {
	Bound result = x_bound(ps, S);
	if (cut.x.valid()) {
		result &= cut.x;
	}
	if (cut.Q_sq.valid()) {
		Bound Q_sq(
			cut.Q_sq.min()/S,
			cut.Q_sq.max()/(cut.Q_sq.max() + sq(ps.Mth + ps.mh) - sq(ps.M)));
		result &= Q_sq;
	}
	if (cut.W_sq.valid()) {
		Real W_sq_min = cut.W_sq.min();
		Real M = ps.M;
		Real m = ps.m;
		Real lambda_S = sq(S) - 4.*sq(M)*sq(m);

		Real L = S - (W_sq_min - sq(M));
		Real denom = 2.*(lambda_S + sq(M)*(W_sq_min - sq(M)));
		Real a = lambda_S - S*(W_sq_min - sq(M));
		Real b = std::sqrt(lambda_S*(sq(L) - 4.*sq(m)*W_sq_min));
		Real x_1 = (a - b)/denom;
		Real x_2 = (a + b)/denom;
		Bound W_sq(x_1, x_2);
		result &= W_sq;
	}
	return result;
}
Bound cut::y_bound(Cut const& cut, Particles const& ps, Real S, Real x) {
	Bound result = y_bound(ps, S, x);
	if (cut.y.valid()) {
		result &= cut.y;
	}
	if (cut.Q_sq.valid()) {
		Bound Q_sq = cut.Q_sq/(S*x);
		result &= Q_sq;
	}
	if (cut.W_sq.valid()) {
		Bound W_sq = (cut.W_sq - sq(ps.M))/((1. - x)*S);
		result &= W_sq;
	}
	return result;
}
Bound cut::z_bound(Cut const& cut, Particles const& ps, Real S, Real x, Real y) {
	Bound result = z_bound(ps, S, x, y);
	if (cut.z.valid()) {
		result &= cut.z;
	}
	return result;
}
Bound cut::ph_t_sq_bound(Cut const& cut, Particles const& ps, Real S, Real x, Real y, Real z) {
	Bound result = ph_t_sq_bound(ps, S, x, y, z);
	if (cut.ph_t_sq.valid()) {
		result &= cut.ph_t_sq;
	}
	return result;
}
Bound cut::tau_bound(CutRad const& cut, Kinematics const& kin) {
	Bound result = tau_bound(kin);
	if (cut.tau.valid()) {
		result &= cut.tau;
	}
	return result;
}
Bound cut::R_bound(CutRad const& cut, Kinematics const& kin, Real tau, Real phi_k) {
	Bound result = R_bound(kin, tau, phi_k);
	if (cut.R.valid()) {
		result &= cut.R;
	}
	if (cut.k_0_bar.valid()) {
		Bound tau_b = tau_bound(kin);
		// Copied from kinematic calculations.
		Real mu = kin.ph_0/kin.M
			+ 1./kin.lambda_Y_sqrt*(
				(2.*tau*sq(kin.M) - kin.S_x)*kin.ph_l/kin.M
				- 2.*kin.M*kin.ph_t*std::cos(kin.phi_h - phi_k)*std::sqrt(
					(tau - tau_b.min())*(tau_b.max() - tau)));
		Bound k_0_bar = (2.*kin.mx*cut.k_0_bar)/(1. + tau - mu);
		result &= k_0_bar;
	}
	return result;
}

// Check whether within kinematic bound.
bool cut::valid(Kinematics const& kin) {
	// TODO: Make sure to check the minimum set of things that let the final
	// state particle be reconstructed.
	if (!(kin.S > 0.)) {
		return false;
	} else if (!(kin.x >= 0. && kin.x <= 1.)) {
		return false;
	} else if (!(kin.y >= 0. && kin.y <= 1.)) {
		return false;
	} else if (!(kin.z >= 0. && kin.z <= 1.)) {
		return false;
	} else if (!(kin.ph_t_sq >= 0.)) {
		return false;
	} else if (!std::isfinite(kin.phi_h)) {
		return false;
	} else if (!std::isfinite(kin.phi)) {
		return false;
	} else if (!(kin.mx >= kin.Mth)) {
		return false;
	} else if (!(kin.mx_sq <= kin.S + sq(kin.M) - sq(kin.mh))) {
		return false;
	} else if (!(kin.q_t >= 0.)) {
		return false;
	} else if (!std::isfinite(kin.q_l)) {
		return false;
	} else if (!(kin.ph_t >= 0.)) {
		return false;
	} else if (!std::isfinite(kin.ph_l)) {
		return false;
	} else if (!(kin.ph_0 >= kin.mh)) {
		return false;
	} else {
		return true;
	}
}
bool cut::valid(KinematicsRad const& kin_rad) {
	// TODO: Instead of checking if things are within bound, check whether the
	// photon can be constructed (which should be equivalent, but more useful).
	if (!valid(kin_rad.project())) {
		return false;
	} else if (!(kin_rad.tau >= kin_rad.tau_min && kin_rad.tau <= kin_rad.tau_max)) {
		return false;
	} else if (!std::isfinite(kin_rad.phi_k)) {
		return false;
	} else if (!(kin_rad.R >= 0 && kin_rad.R <= cut::R_bound(kin_rad.project(), kin_rad.tau, kin_rad.phi_k).max())) {
		return false;
	/*} else if (!(kin_rad.shift_x >= 0. && kin_rad.shift_x <= 1.)) {
		return false;
	} else if (!(kin_rad.shift_z >= 0. && kin_rad.shift_z <= 1.)) {
		return false;
	} else if (!(kin_rad.shift_Q_sq >= 0.)) {
		return false;
	} else if (!(kin_rad.shift_ph_t_sq >= 0.)) {
		return false;*/
	} else {
		return true;
	}
}

bool cut::valid(Cut const& cut, Kinematics const& kin) {
	Real theta_q = std::acos(
		(kin.S*kin.S_x + 2.*sq(kin.M)*kin.Q_sq)
		/(kin.lambda_Y_sqrt*kin.lambda_S_sqrt));
	Real theta_k2 = std::acos(
		(kin.S*kin.X - 2.*sq(kin.M)*kin.Q_sq - 4.*sq(kin.M*kin.m))
		/(kin.lambda_S_sqrt*std::sqrt(sq(kin.X) - 4.*sq(kin.M*kin.m))));
	Real theta_h = std::acos(
		(kin.z*sq(kin.S_x) - 4.*sq(kin.M)*kin.V_m)
		/(kin.lambda_Y_sqrt*std::sqrt(sq(kin.z*kin.S_x) - 4.*sq(kin.M*kin.mh))));
	if (!valid(kin)) {
		return false;
	} else if (cut.x.valid() && !cut.x.contains(kin.x)) {
		return false;
	} else if (cut.y.valid() && !cut.y.contains(kin.y)) {
		return false;
	} else if (cut.z.valid() && !cut.z.contains(kin.z)) {
		return false;
	} else if (cut.ph_t_sq.valid() && !cut.ph_t_sq.contains(kin.ph_t_sq)) {
		return false;
	} else if (cut.phi_h.valid() && !cut.phi_h.contains(kin.phi_h)) {
		return false;
	} else if (cut.phi.valid() && !cut.phi.contains(kin.phi)) {
		return false;
	} else if (cut.Q_sq.valid() && !cut.Q_sq.contains(kin.Q_sq)) {
		return false;
	} else if (cut.t.valid() && !cut.t.contains(kin.t)) {
		return false;
	} else if (cut.W_sq.valid() && !cut.W_sq.contains(kin.W_sq)) {
		return false;
	} else if (cut.r.valid() && !cut.r.contains(kin.V_2/kin.V_1)) {
		return false;
	} else if (cut.mx_sq.valid() && !cut.mx_sq.contains(kin.mx_sq)) {
		return false;
	} else if (cut.q_0.valid() && !cut.q_0.contains(kin.q_0)) {
		return false;
	} else if (cut.k2_0.valid() && !cut.k2_0.contains(kin.k2_0)) {
		return false;
	} else if (cut.ph_0.valid() && !cut.ph_0.contains(kin.ph_0)) {
		return false;
	} else if (cut.theta_q.valid() && !cut.theta_q.contains(theta_q)) {
		return false;
	} else if (cut.theta_k2.valid() && !cut.theta_k2.contains(theta_k2)) {
		return false;
	} else if (cut.theta_h.valid() && !cut.theta_h.contains(theta_h)) {
		return false;
	} else {
		return true;
	}
}

bool cut::valid(CutRad const& cut, KinematicsRad const& kin) {
	Real theta_k = std::acos((kin.S_x - 2.*sq(kin.M)*kin.tau)/kin.lambda_Y_sqrt);
	if (!valid(kin)) {
		return false;
	} else if (cut.tau.valid() && !cut.tau.contains(kin.tau)) {
		return false;
	} else if (cut.phi_k.valid() && !cut.phi_k.contains(kin.phi_k)) {
		return false;
	} else if (cut.R.valid() && !cut.R.contains(kin.R)) {
		return false;
	} else if (cut.k_0_bar.valid() && !cut.k_0_bar.contains(kin.k_0_bar)) {
		return false;
	} else if (cut.k_0.valid() && !cut.k_0.contains(kin.k_0)) {
		return false;
	} else if (cut.theta_k.valid() && !cut.theta_k.contains(theta_k)) {
		return false;
	} else {
		return true;
	}
}

bool cut::valid(Cut const& cut, CutRad const& cut_rad, KinematicsRad const& kin) {
	if (!valid(cut, kin.project())) {
		return false;
	} else if (!valid(cut_rad, kin)) {
		return false;
	} else {
		return true;
	}
}

bool cut::take(
		Particles const& ps, Real S, const Real point[6],
		PhaseSpace* ph_space_out, Real* jac_out) {
	// TODO: A lot of these transformations have been chosen based on how the
	// Born cross-section scales with these variables. However, in certain
	// kinematic regions they might result in very bad integrands--it would be
	// good to investigate these choices further.
	Real jac_x, jac_y, jac_z, jac_ph_t_sq, jac_phi_h, jac_phi;
	Real x = transform(Inverse(X_CUTOFF), point[0], x_bound(ps, S), &jac_x);
	Real y = transform(Inverse(), point[1], y_bound(ps, S, x), &jac_y);
	Real z = transform(Linear(), point[2], z_bound(ps, S, x, y), &jac_z);
	Real ph_t_sq = transform(Decay(PH_T_SQ_DECAY_WIDTH), point[3], ph_t_sq_bound(ps, S, x, y, z), &jac_ph_t_sq);
	Real phi_h = transform(Linear(), point[4], Bound(-PI, PI), &jac_phi_h);
	Real phi = transform(Linear(), point[5], Bound(-PI, PI), &jac_phi);

	if (ph_space_out != nullptr) {
		*ph_space_out = PhaseSpace { x, y, z, ph_t_sq, phi_h, phi };
	}
	if (jac_out != nullptr) {
		*jac_out = jac_x * jac_y * jac_z * jac_ph_t_sq * jac_phi_h * jac_phi;
	}
	// TODO: Should we check for validity explicitly here?
	return true;
}

bool cut::take(
		Particles const& ps, Real S, const Real point[6],
		Kinematics* kin_out, Real* jac_out) {
	PhaseSpace ph_space;
	if (!cut::take(ps, S, point, &ph_space, jac_out)) {
		return false;
	}
	if (kin_out != nullptr) {
		*kin_out = Kinematics(ps, S, ph_space);
	}
	return true;
}

bool cut::take(
		Cut const& cut,
		Particles const& ps, Real S, const Real point[6],
		Kinematics* kin_out, Real* jac_out) {
	Real jac_x, jac_y, jac_z, jac_ph_t_sq, jac_phi_h, jac_phi;
	Real x = transform(Inverse(X_CUTOFF), point[0], x_bound(cut, ps, S), &jac_x);
	Real y = transform(Inverse(), point[1], y_bound(cut, ps, S, x), &jac_y);
	Real z = transform(Linear(), point[2], z_bound(cut, ps, S, x, y), &jac_z);
	Real ph_t_sq = transform(Decay(PH_T_SQ_DECAY_WIDTH), point[3], ph_t_sq_bound(cut, ps, S, x, y, z), &jac_ph_t_sq);
	Real phi_h = transform(Linear(), point[4], cut.phi_h.valid() ? cut.phi_h : Bound(-PI, PI), &jac_phi_h);
	Real phi = transform(Linear(), point[5], cut.phi.valid() ? cut.phi : Bound(-PI, PI), &jac_phi);

	PhaseSpace ph_space { x, y, z, ph_t_sq, phi_h, phi };
	Kinematics kin(ps, S, ph_space);
	if (kin_out != nullptr) {
		*kin_out = kin;
	}
	if (jac_out != nullptr) {
		*jac_out = jac_x * jac_y * jac_z * jac_ph_t_sq * jac_phi_h * jac_phi;
	}
	return valid(cut, kin);
}

bool cut::take(
		Particles const& ps, Real S, const Real point[9],
		PhaseSpaceRad* ph_space_out, Real* jac_out) {
	Kinematics kin;
	Real jac;
	if (!take(ps, S, point, &kin, &jac)) {
		return false;
	}
	if (!take(kin, point + 6, ph_space_out, jac_out)) {
		return false;
	}
	if (jac_out != nullptr) {
		*jac_out *= jac;
	}
	return true;
}

bool cut::take(
		Particles const& ps, Real S, const Real point[9],
		KinematicsRad* kin_out, Real* jac_out) {
	Kinematics kin;
	Real jac;
	if (!take(ps, S, point, &kin, &jac)) {
		return false;
	}
	if (!take(kin, point + 6, kin_out, jac_out)) {
		return false;
	}
	if (jac_out != nullptr) {
		*jac_out *= jac;
	}
	return true;
}

bool cut::take(
		Cut const& cut, CutRad const& cut_rad,
		Particles const& ps, Real S, const Real point[9],
		KinematicsRad* kin_out, Real* jac_out) {
	Kinematics kin;
	Real jac;
	if (!take(cut, ps, S, point, &kin, &jac)) {
		return false;
	}
	if (!take(cut_rad, kin, point + 6, kin_out, jac_out)) {
		return false;
	}
	if (jac_out != nullptr) {
		*jac_out *= jac;
	}
	return true;
}

bool cut::take(
		Kinematics const& kin, const Real point[3],
		PhaseSpaceRad* ph_space_out, Real* jac_out) {
	// Calculate peaks in radiative distribution.
	Real tau_p1 = kin.Q_sq / kin.X;
	Real tau_pr = -(1. - kin.y);
	Real tau_lim1 = SQRT_2 * kin.lambda_1_sqrt * kin.m / sq(kin.X);
	Real tau_limr = sq(1. - kin.y);
	Real phi_k_lim = kin.lambda_Y_sqrt / kin.lambda_1_sqrt / SQRT_2 * kin.m;
	// Estimate of the point at which `R` starts behaving like inverse `tau`.
	// TODO: The factor of 128 is arbitrary right now, try to tune it.
	Real R_trans = (kin.mx_sq - sq(kin.Mth))
		/ (1. + kin.S_x / sq(kin.M) - kin.ph_l / kin.M) / 128.;

	Real jac_tau, jac_phi_k, jac_R;
	Real tau = transform(
		Sigmoid2(tau_p1, tau_lim1, tau_pr, tau_limr),
		point[0], tau_bound(kin), &jac_tau);
	Real phi_k = transform(
		Sigmoid(0., phi_k_lim),
		point[1], Bound(-PI, PI), &jac_phi_k);
	Real R = transform(
		Log(R_trans),
		point[2], R_bound(kin, tau, phi_k), &jac_R);

	if (ph_space_out != nullptr) {
		*ph_space_out = {
			kin.x, kin.y, kin.z,
			kin.ph_t_sq, kin.phi_h, kin.phi,
			tau, phi_k, R,
		};
	}
	if (jac_out != nullptr) {
		*jac_out = jac_tau * jac_phi_k * jac_R;
	}
	// TODO: Should we check explicitly for validity here?
	return true;
}

bool cut::take(
		Kinematics const& kin, const Real point[3],
		KinematicsRad* kin_out, Real* jac_out) {
	PhaseSpaceRad ph_space;
	if (!cut::take(kin, point, &ph_space, jac_out)) {
		return false;
	}
	if (kin_out != nullptr) {
		*kin_out = KinematicsRad(kin, ph_space.tau, ph_space.phi_k, ph_space.R);
	}
	return true;
}

bool cut::take(
		CutRad const& cut,
		Kinematics const& kin, const Real point[3],
		KinematicsRad* kin_out, Real* jac_out) {
	Real tau_p1 = kin.Q_sq / kin.X;
	Real tau_pr = -(1. - kin.y);
	Real tau_lim1 = SQRT_2 * kin.lambda_1_sqrt * kin.m / sq(kin.X);
	Real tau_limr = sq(1. - kin.y);
	Real phi_k_lim = kin.lambda_Y_sqrt / kin.lambda_1_sqrt / SQRT_2 * kin.m;
	Real R_trans = (kin.mx_sq - sq(kin.Mth))
		/ (1. + kin.S_x / sq(kin.M) - kin.ph_l / kin.M) / 128.;

	Real jac_tau, jac_phi_k, jac_R;
	Real tau = transform(
		Sigmoid2(tau_p1, tau_lim1, tau_pr, tau_limr),
		point[0], tau_bound(cut, kin), &jac_tau);
	Real phi_k = transform(
		Sigmoid(0., phi_k_lim),
		point[1], cut.phi_k.valid() ? cut.phi_k : Bound(-PI, PI), &jac_phi_k);
	Real R = transform(
		Log(R_trans),
		point[2], R_bound(cut, kin, tau, phi_k), &jac_R);

	KinematicsRad kin_rad(kin, tau, phi_k, R);
	if (kin_out != nullptr) {
		*kin_out = kin_rad;
	}
	if (jac_out != nullptr) {
		*jac_out = jac_tau * jac_phi_k * jac_R;
	}
	return valid(cut, kin_rad);
}

