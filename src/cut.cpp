#include "sidis/cut.hpp"

#include <algorithm>
#include <cmath>

#include "sidis/constant.hpp"
#include "sidis/kinematics.hpp"
#include "sidis/extra/math.hpp"

using namespace sidis;
using namespace sidis::constant;
using namespace sidis::cut;
using namespace sidis::kin;
using namespace sidis::math;

Cut::Cut() :
	x(Bounds::UNIT),
	y(Bounds::UNIT),
	z(Bounds::UNIT),
	ph_t_sq(Bounds::POSITIVE),
	phi_h(Bounds(-PI, PI)),
	phi(Bounds(-PI, PI)) { }

CutRad::CutRad() :
	tau(Bounds::FULL),
	phi_k(Bounds(-PI, PI)),
	k_0_bar(Bounds::POSITIVE) { }

// Bounds of kinematic variables.
// TODO: Some of the calculations in this section are redundant with earlier
// kinematic calculations. This should be refactored to avoid that later.
Real cut::S_min(Particles ps) {
	return sq(ps.Mth + ps.mh) + 2.*ps.m*(ps.Mth + ps.mh) - sq(ps.M);
}
Bounds cut::x_bounds(Particles ps, Real S) {
	Real M = ps.M;
	Real m = ps.m;
	Real mh = ps.mh;
	Real Mth = ps.Mth;
	Real lambda_S = sq(S) - 4.*sq(M)*sq(m);
	Real L = S - (sq(Mth + mh) - sq(M));

	// Combined bounds `mx >= Mth` and `sq(q_t) >= 0`.
	Real denom = 2.*(lambda_S + sq(M)*(sq(Mth + mh) - sq(M)));
	Real a = lambda_S - S*(sq(Mth + mh) - sq(M));
	Real b = std::sqrt(lambda_S*(L - 2.*m*(Mth + mh))*(L + 2.*m*(Mth + mh)));
	Real x_1 = (a - b)/denom;
	Real x_2 = (a + b)/denom;
	Bounds kin_b(x_1, x_2);

	return Bounds::UNIT & kin_b;
}
Bounds cut::y_bounds(Particles ps, Real S, Real x) {
	Real M = ps.M;
	Real m = ps.m;
	Real mh = ps.mh;
	Real Mth = ps.Mth;
	Real lambda_S = sq(S) - 4.*sq(M)*sq(m);

	// Bound `mx >= Mth`.
	Real px_threshold = (sq(Mth + mh) - sq(M))/((1. - x)*S);
	// Bound `sq(q_t) >= 0`.
	Real q_threshold = (x*lambda_S)/(S*(sq(M*x) + S*x + sq(m)));
	Bounds kin_b(px_threshold, q_threshold);

	return Bounds::UNIT & kin_b;
}
Bounds cut::z_bounds(Particles ps, Real S, Real x, Real y) {
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
	Bounds kin_b(px_threshold_min, px_threshold_max);

	return Bounds::UNIT & kin_b;
}
Bounds cut::ph_t_sq_bounds(Particles ps, Real S, Real x, Real y, Real z) {
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
	Bounds kin_b(0., px_threshold);

	return Bounds::POSITIVE & kin_b;
}
Bounds cut::tau_bounds(Kinematics kin) {
	// Equation [1.44].
	Real min = (kin.S_x - kin.lambda_Y_sqrt)/(2.*sq(kin.M));
	Real max = (kin.S_x + kin.lambda_Y_sqrt)/(2.*sq(kin.M));
	Bounds kin_b(min, max);

	return Bounds::FULL & kin_b;
}
Bounds cut::R_bounds(Kinematics kin, Real tau, Real phi_k) {
	Bounds tau_b = tau_bounds(kin);
	// Copied from kinematic calculations.
	Real mu = kin.ph_0/kin.M
		+ 1./kin.lambda_Y_sqrt*(
			(2.*tau*sq(kin.M) - kin.S_x)*kin.ph_l/kin.M
			- 2.*kin.M*kin.ph_t*std::cos(kin.phi_h - phi_k)*std::sqrt(
				(tau - tau_b.min())*(tau_b.max() - tau)));
	// Equation [1.44].
	Real R_max = 1./(1. + tau - mu)*(kin.mx_sq - sq(kin.Mth));
	Bounds kin_b(0., R_max);

	return Bounds::POSITIVE & kin_b;
}

Bounds cut::x_bounds(Cut cut, Particles ps, Real S) {
	return x_bounds(ps, S) & cut.x;
}
Bounds cut::y_bounds(Cut cut, Particles ps, Real S, Real x) {
	return y_bounds(ps, S, x) & cut.y;
}
Bounds cut::z_bounds(Cut cut, Particles ps, Real S, Real x, Real y) {
	return z_bounds(ps, S, x, y) & cut.z;
}
Bounds cut::ph_t_sq_bounds(Cut cut, Particles ps, Real S, Real x, Real y, Real z) {
	return ph_t_sq_bounds(ps, S, x, y, z) & cut.ph_t_sq;
}
Bounds cut::tau_bounds(CutRad cut, Kinematics kin) {
	return tau_bounds(kin) & cut.tau;
}
Bounds cut::R_bounds(CutRad cut, Kinematics kin, Real tau, Real phi_k) {
	Bounds tau_b = tau_bounds(kin);
	// Copied from kinematic calculations.
	Real mu = kin.ph_0/kin.M
		+ 1./kin.lambda_Y_sqrt*(
			(2.*tau*sq(kin.M) - kin.S_x)*kin.ph_l/kin.M
			- 2.*kin.M*kin.ph_t*std::cos(kin.phi_h - phi_k)*std::sqrt(
				(tau - tau_b.min())*(tau_b.max() - tau)));
	Bounds k_0_bar_b = 1./(1. + tau - mu)*(2.*kin.mx*cut.k_0_bar - sq(kin.Mth));
	return R_bounds(kin, tau, phi_k) & k_0_bar_b;
}

// Check whether within kinematic bounds.
bool cut::valid(Kinematics kin) {
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
bool cut::valid(KinematicsRad kin_rad) {
	// TODO: Instead of checking if things are within bounds, check whether the
	// photon can be constructed (which should be equivalent, but more useful).
	if (!valid(kin_rad.project())) {
		return false;
	} else if (!(kin_rad.tau >= kin_rad.tau_min && kin_rad.tau <= kin_rad.tau_max)) {
		return false;
	} else if (!std::isfinite(kin_rad.phi_k)) {
		return false;
	} else if (!(kin_rad.R >= 0. && kin_rad.R <= kin_rad.R_max)) {
		return false;
	} else {
		return true;
	}
}

bool cut::valid(Cut cut, Kinematics kin) {
	if (!valid(kin)) {
		return false;
	} else if (!cut.x.contains(kin.x)) {
		return false;
	} else if (!cut.y.contains(kin.y)) {
		return false;
	} else if (!cut.z.contains(kin.z)) {
		return false;
	} else if (!cut.ph_t_sq.contains(kin.ph_t_sq)) {
		return false;
	} else if (!cut.phi_h.contains(kin.phi_h)) {
		return false;
	} else if (!cut.phi.contains(kin.phi)) {
		return false;
	} else {
		return true;
	}
}

bool cut::valid(CutRad cut, KinematicsRad kin) {
	if (!valid(kin)) {
		return false;
	} else if (!cut.tau.contains(kin.tau)) {
		return false;
	} else if (!cut.phi_k.contains(kin.phi_k)) {
		return false;
	} else if (!cut.k_0_bar.contains(kin.k_0_bar)) {
		return false;
	} else {
		return true;
	}
}

bool cut::valid(Cut cut, CutRad cut_rad, KinematicsRad kin) {
	if (!valid(cut, kin.project())) {
		return false;
	} else if (!valid(cut_rad, kin)) {
		return false;
	} else {
		return true;
	}
}

bool cut::take(
		Particles ps, Real S, const Real point[6],
		PhaseSpace* ph_space_out, Real* jacobian_out) {
	Bounds x_b = x_bounds(ps, S);
	Real x = x_b.lerp(point[0]);
	Bounds y_b = y_bounds(ps, S, x);
	Real y = y_b.lerp(point[1]);
	Bounds z_b = z_bounds(ps, S, x, y);
	Real z = z_b.lerp(point[2]);
	Bounds ph_t_sq_b = ph_t_sq_bounds(ps, S, x, y, z);
	Real ph_t_sq = ph_t_sq_b.lerp(point[3]);
	Bounds phi_h_b = Bounds(-PI, PI);
	Real phi_h = phi_h_b.lerp(point[4]);
	Bounds phi_b = Bounds(-PI, PI);
	Real phi = phi_b.lerp(point[5]);
	if (ph_space_out != nullptr) {
		*ph_space_out = PhaseSpace { x, y, z, ph_t_sq, phi_h, phi };
	}
	if (jacobian_out != nullptr) {
		*jacobian_out = x_b.size() * y_b.size() * z_b.size()
			* ph_t_sq_b.size() * phi_h_b.size() * phi_b.size();
	}
	// TODO: Should we check for validity explicitly here?
	return true;
}

bool cut::take(
		Particles ps, Real S, const Real point[6],
		Kinematics* kin_out, Real* jacobian_out) {
	PhaseSpace ph_space;
	if (!cut::take(ps, S, point, &ph_space, jacobian_out)) {
		return false;
	}
	if (kin_out != nullptr) {
		*kin_out = Kinematics(ps, S, ph_space);
	}
	return true;
}

bool cut::take(
		Cut cut, Particles ps, Real S, const Real point[6],
		Kinematics* kin_out, Real* jacobian_out) {
	Bounds x_b = x_bounds(cut, ps, S);
	Real x = x_b.lerp(point[0]);
	Bounds y_b = y_bounds(cut, ps, S, x);
	Real y = y_b.lerp(point[1]);
	Bounds z_b = z_bounds(cut, ps, S, x, y);
	Real z = z_b.lerp(point[2]);
	Bounds ph_t_sq_b = ph_t_sq_bounds(cut, ps, S, x, y, z);
	Real ph_t_sq = ph_t_sq_b.lerp(point[3]);
	Bounds phi_h_b = cut.phi_h;
	Real phi_h = phi_h_b.lerp(point[4]);
	Bounds phi_b = cut.phi;
	Real phi = phi_b.lerp(point[5]);
	PhaseSpace ph_space { x, y, z, ph_t_sq, phi_h, phi };
	Kinematics kin(ps, S, ph_space);
	if (kin_out != nullptr) {
		*kin_out = kin;
	}
	if (jacobian_out != nullptr) {
		*jacobian_out = x_b.size() * y_b.size() * z_b.size()
			* ph_t_sq_b.size() * phi_h_b.size() * phi_b.size();
	}
	return valid(cut, kin);
}

bool cut::take(
		Particles ps, Real S, const Real point[9],
		PhaseSpaceRad* ph_space_out, Real* jacobian_out) {
	Kinematics kin;
	Real jacobian;
	if (!take(ps, S, point, &kin, &jacobian)) {
		return false;
	}
	if (!take(kin, point + 6, ph_space_out, jacobian_out)) {
		return false;
	}
	if (jacobian_out != nullptr) {
		*jacobian_out *= jacobian;
	}
	return true;
}

bool cut::take(
		Particles ps, Real S, const Real point[9],
		KinematicsRad* kin_out, Real* jacobian_out) {
	Kinematics kin;
	Real jacobian;
	if (!take(ps, S, point, &kin, &jacobian)) {
		return false;
	}
	if (!take(kin, point + 6, kin_out, jacobian_out)) {
		return false;
	}
	if (jacobian_out != nullptr) {
		*jacobian_out *= jacobian;
	}
	return true;
}

bool cut::take(
		Cut cut, CutRad cut_rad, Particles ps, Real S, const Real point[9],
		KinematicsRad* kin_out, Real* jacobian_out) {
	Kinematics kin;
	Real jacobian;
	if (!take(cut, ps, S, point, &kin, &jacobian)) {
		return false;
	}
	if (!take(cut_rad, kin, point + 6, kin_out, jacobian_out)) {
		return false;
	}
	if (jacobian_out != nullptr) {
		*jacobian_out *= jacobian;
	}
	return true;
}

bool cut::take(
		Kinematics kin, const Real point[3],
		PhaseSpaceRad* ph_space_out, Real* jacobian_out) {
	Bounds tau_b = tau_bounds(kin);
	Real tau = tau_b.lerp(point[6]);
	Bounds phi_k_b = Bounds(-PI, PI);
	Real phi_k = phi_k_b.lerp(point[7]);
	Bounds R_b = R_bounds(kin, tau, phi_k);
	Real R = R_b.lerp(point[8]);
	if (ph_space_out != nullptr) {
		*ph_space_out = {
			kin.x, kin.y, kin.z,
			kin.ph_t_sq, kin.phi_h,
			kin.phi, tau, phi_k, R,
		};
	}
	if (jacobian_out != nullptr) {
		*jacobian_out = tau_b.size() * phi_k_b.size() * R_b.size();
	}
	// TODO: Should we check explicitly for validity here?
	return true;
}

bool cut::take(
		Kinematics kin, const Real point[3],
		KinematicsRad* kin_out, Real* jacobian_out) {
	PhaseSpaceRad ph_space;
	if (!cut::take(kin, point, &ph_space, jacobian_out)) {
		return false;
	}
	if (kin_out != nullptr) {
		*kin_out = KinematicsRad(kin, ph_space.tau, ph_space.phi_k, ph_space.R);
	}
	return true;
}

bool cut::take(
		CutRad cut, Kinematics kin, const Real point[3],
		KinematicsRad* kin_out, Real* jacobian_out) {
	Bounds tau_b = tau_bounds(cut, kin);
	Real tau = tau_b.lerp(point[6]);
	Bounds phi_k_b = cut.phi_k;
	Real phi_k = phi_k_b.lerp(point[7]);
	Bounds R_b = R_bounds(cut, kin, tau, phi_k);
	Real R = R_b.lerp(point[8]);
	KinematicsRad kin_rad(kin, tau, phi_k, R);
	if (kin_out != nullptr) {
		*kin_out = kin_rad;
	}
	if (jacobian_out != nullptr) {
		*jacobian_out = tau_b.size() * phi_k_b.size() * R_b.size();
	}
	return valid(cut, kin_rad);
}

