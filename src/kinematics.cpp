#include "sidis/kinematics.hpp"

#include <cmath>

#include "sidis/constant.hpp"
#include "sidis/frame.hpp"
#include "sidis/extra/bounds.hpp"
#include "sidis/extra/math.hpp"
#include "sidis/extra/transform.hpp"

using namespace sidis;
using namespace sidis::constant;
using namespace sidis::frame;
using namespace sidis::kin;
using namespace sidis::math;

Kinematics::Kinematics(
		Initial init,
		PhaseSpace ph_space,
		Hadron hadron,
		Real M_th) {
	x = ph_space.x;
	y = ph_space.y;
	z = ph_space.z;
	ph_t_sq = ph_space.ph_t_sq;
	phi_h = ph_space.phi_h;
	phi = ph_space.phi;
	phi_q = std::fmod(PI - ph_space.phi, 2.*PI);
	if (phi_q > PI) {
		phi_q -= 2.*PI;
	} else if (phi_q < -PI) {
		phi_q += 2.*PI;
	}

	cos_phi_h = std::cos(phi_h);
	sin_phi_h = std::sin(phi_h);
	cos_phi = std::cos(phi);
	sin_phi = std::sin(phi);
	cos_phi_q = -cos_phi;
	sin_phi_q = sin_phi;

	target = init.target;
	beam = init.beam;
	this->hadron = hadron;

	S = 2.*dot(init.p, init.k1);
	M = mass(target);
	m = mass(beam);
	mh = mass(hadron);
	this->M_th = M_th;

	// Equation [1.3].
	Q_sq = S*x*y;
	Q = std::sqrt(Q_sq);
	X = S*(1. - y);
	S_x = S*y;
	S_p = S*(2. - y);
	lambda_S = sq(S) - 4.*sq(M)*sq(m);
	lambda_Y = sq(S_x) + 4.*sq(M)*Q_sq;
	lambda_1 = Q_sq*(S*X - sq(M)*Q_sq) - sq(m)*lambda_Y;
	lambda_S_sqrt = std::sqrt(lambda_S);
	lambda_Y_sqrt = std::sqrt(lambda_Y);
	lambda_1_sqrt = std::sqrt(lambda_1);

	// Equation [1.4]. The equations have been re-arranged in terms of
	// `ph_t_sq`.
	ph_0 = (z*S_x)/(2.*M);
	ph_t = std::sqrt(ph_t_sq);
	Real ph_ratio_sq = ph_t_sq/sq(ph_0) + sq(mh/ph_0);
	ph_l = ph_0*std::sqrt(1. - ph_ratio_sq);

	// Virtual photon 4-momentum components.
	q_0 = S_x/(2.*M);
	// Equation [1.4].
	q_t = lambda_1_sqrt/lambda_S_sqrt;
	q_l = (2.*sq(M)*Q_sq + S*S_x)/(2.*M*lambda_S_sqrt);

	// Scattered lepton 4-momentum components.
	k2_0 = X/(2.*M);
	k2_l = (S*X - 2.*sq(M)*Q_sq - 4.*sq(M)*sq(m))/(2.*M*lambda_S_sqrt);
	k2_t = lambda_1_sqrt/lambda_S_sqrt;
	k1_t = lambda_1_sqrt/lambda_Y_sqrt;

	// Equation [1.5].
	V_1 = ph_0*S/M - (ph_l*(S*S_x + 2.*sq(M)*Q_sq))/(M*lambda_Y_sqrt)
		- 2.*ph_t*k1_t*cos_phi_h;
	V_2 = ph_0*X/M - (ph_l*(X*S_x - 2.*sq(M)*Q_sq))/(M*lambda_Y_sqrt)
		- 2.*ph_t*k1_t*cos_phi_h;
	V_p = 0.5*(V_1 + V_2);

	// In the low `ph_t` case (where the cross-section is the highest), the
	// computation for `V_m` has a catastrophic cancellation between the terms
	// `2 M ph_l √λ_Y - z S_x²`. So, it's better to compute it in the following
	// way:
	Real lambda_Y_ratio = (4.*sq(M)*Q_sq)/sq(S_x);
	V_m = -(ph_0*S_x)/(2.*M)*sqrt1p_1m(
		lambda_Y_ratio - ph_ratio_sq - lambda_Y_ratio*ph_ratio_sq);

	t = sq(mh) - Q_sq - 2.*V_m;
	mx_sq = sq(M) + t + (1. - z)*S_x;
	mx = std::sqrt(mx_sq);

	// Paragraph below equation [1.14].
	lambda_2 = sq(V_m) + sq(mh)*Q_sq;
	lambda_3 = V_m + z*Q_sq;
	lambda_2_sqrt = std::sqrt(lambda_2);
	lambda_3_sqrt = std::sqrt(lambda_3);

	// Equation [1.6]. `vol_phi_h` is defined as `dot(epsilon_perp, ph)`.
	vol_phi_h = -0.5*ph_t*lambda_1_sqrt*sin_phi_h;

	// Equation [1.18]. Modified slightly to include a `Q^4` denominator.
	C_1 = (4.*M*ph_l*(Q_sq + 2.*x*sq(M)))/sq(sq(Q_sq));
}

KinematicsRad::KinematicsRad(Kinematics kin, Real tau, Real phi_k, Real R) :
		target(kin.target),
		beam(kin.beam),
		hadron(kin.hadron),
		S(kin.S),
		M(kin.M),
		m(kin.m),
		mh(kin.mh),
		M_th(kin.M_th),
		x(kin.x),
		y(kin.y),
		z(kin.z),
		ph_t_sq(kin.ph_t_sq),
		phi_h(kin.phi_h),
		phi(kin.phi),
		phi_q(kin.phi_q),
		cos_phi_h(kin.cos_phi_h),
		sin_phi_h(kin.sin_phi_h),
		cos_phi(kin.cos_phi),
		sin_phi(kin.sin_phi),
		cos_phi_q(kin.cos_phi_q),
		sin_phi_q(kin.sin_phi_q),
		Q_sq(kin.Q_sq),
		Q(kin.Q),
		t(kin.t),
		X(kin.X),
		S_x(kin.S_x),
		S_p(kin.S_p),
		V_1(kin.V_1),
		V_2(kin.V_2),
		V_m(kin.V_m),
		V_p(kin.V_p),
		lambda_S(kin.lambda_S),
		lambda_Y(kin.lambda_Y),
		lambda_1(kin.lambda_1),
		lambda_2(kin.lambda_2),
		lambda_3(kin.lambda_3),
		lambda_S_sqrt(kin.lambda_S_sqrt),
		lambda_Y_sqrt(kin.lambda_Y_sqrt),
		lambda_1_sqrt(kin.lambda_1_sqrt),
		lambda_2_sqrt(kin.lambda_2_sqrt),
		lambda_3_sqrt(kin.lambda_3_sqrt),
		ph_0(kin.ph_0),
		ph_t(kin.ph_t),
		ph_l(kin.ph_l),
		q_0(kin.q_0),
		q_t(kin.q_t),
		q_l(kin.q_l),
		k2_0(kin.k2_0),
		k2_t(kin.k2_t),
		k2_l(kin.k2_l),
		k1_t(kin.k1_t),
		mx_sq(kin.mx_sq),
		mx(kin.mx),
		vol_phi_h(kin.vol_phi_h),
		C_1(kin.C_1),
		tau(tau),
		phi_k(phi_k),
		R(R) {
	cos_phi_k = std::cos(phi_k);
	sin_phi_k = std::sin(phi_k);

	// Equation [1.44].
	tau_min = (S_x - lambda_Y_sqrt)/(2.*sq(M));
	tau_max = (S_x + lambda_Y_sqrt)/(2.*sq(M));

	// TODO: Fill in equation number from derivations.
	lambda_H = sq(z*S_x) - 4.*sq(M)*sq(mh);
	// TODO: Fill in equation number from derivations.
	lambda_V = z*sq(S_x) - 4.*sq(M)*V_m;
	lambda_RY = R*(S_x - 2.*sq(M)*tau);
	lambda_RV = (2.*M)/lambda_Y_sqrt*(
		2.*sq(M)*R*ph_t*std::sqrt((tau - tau_min)*(tau_max - tau))
			*std::cos(phi_h - phi_k)
		+ lambda_RY*ph_l);

	// Equation [1.B3]. The equation has been modified to account for our sign
	// conventions on the angles `phi_h` and `phi_k`.
	mu = ph_0/M
		+ 1./lambda_Y_sqrt*(
			(2.*tau*sq(M) - S_x)*ph_l/M
			- 2.*M*ph_t*std::cos(phi_h - phi_k)*std::sqrt(
				(tau - tau_min)*(tau_max - tau)));
	// Alternative definition of `mu`. Undefined at `R=0`.
	//   mu = (z*R*S_x - lambda_RV)/(2.*sq(M)*R);

	// Equation [1.44].
	R_max = (mx_sq - sq(M_th))/(1. + tau - mu);

	// Equation [1.B4].
	lambda_z = (tau_max - tau)*(tau - tau_min)*lambda_1;
	lambda_z_sqrt = std::sqrt(lambda_z);
	z_1 = 1./lambda_Y*(
		Q_sq*S_p
		+ tau*(S*S_x + 2.*sq(M)*Q_sq)
		- 2.*M*lambda_z_sqrt*cos_phi_k);
	z_2 = 1./lambda_Y*(
		Q_sq*S_p
		+ tau*(X*S_x - 2.*sq(M)*Q_sq)
		- 2.*M*lambda_z_sqrt*cos_phi_k);

	// Real photon 4-momentum components.
	k_0 = R/(2.*M);
	k_t = (M*R*lambda_z_sqrt)/(lambda_1_sqrt*lambda_Y_sqrt);
	k_l = lambda_RY/(2.*M*lambda_Y_sqrt);

	// Equation [1.B5].
	F_22 = 1./sq(z_2);
	F_21 = 1./sq(z_1);
	F_2p = F_22 + F_21;
	F_2m = F_22 - F_21;
	F_d = 1./(z_1*z_2);
	F_1p = 1./z_1 + 1./z_2;
	F_IR = sq(m)*F_2p - (Q_sq + 2.*sq(m))*F_d;

	// Equation [1.30]. Modified to remove the factor of `R`.
	// TODO: Why does this equation require a negative sign compared to what is
	// given in [1]? One possibility is that we use an opposite sign convention
	// than what is used in [1] for `phi_k`.
	vol_phi_k_R = -0.5*sin_phi_k*M*lambda_z_sqrt/lambda_Y_sqrt;
	// Equation [1.A9].
	vol_phi_hk = 1./(2.*lambda_1)*(
		R*vol_phi_h*(z_1*lambda_Y - Q_sq*S_p - tau*(S*S_x + 2.*sq(M)*Q_sq))
		+ R*vol_phi_k_R*(S_x*(z*Q_sq*S_p - S*V_2 + X*V_1) - 4.*V_p*sq(M)*Q_sq));

	shift_Q_sq = Q_sq + R*tau;
	shift_Q = std::sqrt(shift_Q_sq);
	shift_S_x = S_x - R;
	shift_V_m = V_m + (lambda_RV - z*R*S_x)/(4.*sq(M));

	shift_x = shift_Q_sq/shift_S_x;
	shift_y = y - R/S;
	shift_z = (2.*M*ph_0)/shift_S_x;
	// TODO: Fill in equation number in derivations.
	shift_t = t - R*tau + (z*R*S_x - lambda_RV)/(2.*sq(M));

	shift_lambda_Y = lambda_Y + sq(R) - 2.*lambda_RY;
	shift_lambda_1 = lambda_1 + 1./(4.*sq(M))*(
		(sq(R) - 2.*lambda_RY)*lambda_S
		+ R*(S - 2.*sq(M)*z_1)*(
			2.*S*S_x + 4.*sq(M)*Q_sq
			- R*(S - 2.*sq(M)*z_1)));
	shift_lambda_2 = sq(shift_V_m) + sq(mh)*shift_Q_sq;
	shift_lambda_3 = shift_V_m + shift_z*shift_Q_sq;
	shift_lambda_Y_sqrt = std::sqrt(shift_lambda_Y);
	shift_lambda_1_sqrt = std::sqrt(shift_lambda_1);
	shift_lambda_2_sqrt = std::sqrt(shift_lambda_2);
	shift_lambda_3_sqrt = std::sqrt(shift_lambda_3);

	shift_ph_t_sq = ph_t_sq + 1./(shift_lambda_Y)*(
		+ (sq(R) - 2.*lambda_RY)*sq(ph_l)
		+ (lambda_Y_sqrt*lambda_RV*ph_l)/M
		- sq(lambda_RV)/(4.*sq(M)));
	shift_ph_t = std::sqrt(shift_ph_t_sq);
	shift_ph_l = 1./shift_lambda_Y_sqrt*(lambda_Y_sqrt*ph_l - lambda_RV/(2.*M));
	shift_q_0 = shift_S_x/(2.*M);
	// TODO: Find a more accurate method for calculating `shift_q_t`.
	shift_q_t = std::sqrt(
		sq(q_t)
		+ R/(4.*sq(M))*(R - 2.*(S_x - 2.*sq(M)*tau))
		- R/(4.*sq(M)*lambda_S)*(S - 2.*sq(M)*z_1)*(
			R*(S - 2.*sq(M)*z_1) - 2.*(S*S_x + 2.*sq(M)*Q_sq)));
	shift_q_l = q_l - R/(2.*M*lambda_S_sqrt)*(S - 2.*sq(M)*z_1);
	shift_k1_t = shift_lambda_1_sqrt/shift_lambda_Y_sqrt;
	shift_mx_sq = mx_sq - R*(1. + tau) + (z*R*S_x - lambda_RV)/(2.*sq(M));
	shift_mx = std::sqrt(shift_mx_sq);

	// Equation [1.A9].
	shift_vol_phi_h = vol_phi_h + 1./(2.*lambda_1)*(
		R*vol_phi_h*(
			tau*lambda_S
			+ 2.*sq(m)*S_x
			+ Q_sq*S
			- z_1*(S*S_x + 2.*sq(M)*Q_sq))
		+ R*vol_phi_k_R*(
			2.*sq(m)*(4.*V_m*sq(M) - z*sq(S_x))
			+ S*(S*V_2 - X*V_1 - z*Q_sq*S_x)
			+ 2.*V_1*sq(M)*Q_sq));

	shift_C_1 = (4.*M*shift_ph_l*(shift_Q_sq + 2.*shift_x*sq(M)))/sq(sq(shift_Q_sq));

	// TODO: Fill in equation number from derivations.
	shift_sin_phi_h = -2.*shift_vol_phi_h/(shift_ph_t*shift_q_t*lambda_S_sqrt);
	shift_cos_phi_h = 1./(4.*sq(M)*shift_ph_t*shift_q_t*shift_lambda_Y_sqrt*lambda_S_sqrt)*(
		shift_lambda_Y*(z*S*S_x - 2.*sq(M)*V_1)
		- (lambda_V - lambda_RV)*(S*shift_S_x + 2.*sq(M)*Q_sq + 2.*sq(M)*z_1*R));
	shift_phi_h = std::atan2(shift_sin_phi_h, shift_cos_phi_h);

	// TODO: Verify the correctness of these expressions (considering the
	// opposite signs of the `phi` and `phi_q` angles).
	shift_sin_phi_q = 1./shift_q_t*(
		q_t*sin_phi_q
		- (2.*M)/lambda_Y_sqrt*(
			k_t*(-q_l*sin_phi*cos_phi_k + lambda_Y_sqrt/(2.*M)*cos_phi*sin_phi_k)
			+ k_l*q_t*sin_phi));
	shift_cos_phi_q = 1./shift_q_t*(
		q_t*cos_phi_q
		- (2.*M)/lambda_Y_sqrt*(
			k_t*(q_l*cos_phi*cos_phi_k + lambda_Y_sqrt/(2.*M)*sin_phi*sin_phi_k)
			- k_l*q_t*cos_phi));
	shift_phi_q = std::atan2(shift_sin_phi_q, shift_cos_phi_q);
}

Kinematics KinematicsRad::project() const {
	Kinematics kin;
	kin.target = target;
	kin.beam = beam;
	kin.hadron = hadron;

	kin.S = S;
	kin.M = M;
	kin.m = m;
	kin.mh = mh;
	kin.M_th = M_th;

	kin.x = x;
	kin.y = y;
	kin.z = z;
	kin.ph_t_sq = ph_t_sq;
	kin.phi_h = phi_h;
	kin.phi = phi;

	kin.phi_q = phi_q;

	kin.cos_phi_h = cos_phi_h;
	kin.sin_phi_h = sin_phi_h;
	kin.cos_phi = cos_phi;
	kin.sin_phi = sin_phi;
	kin.cos_phi_q = cos_phi_q;
	kin.sin_phi_q = sin_phi_q;

	kin.Q_sq = Q_sq;
	kin.Q = Q;
	kin.t = t;
	kin.X = X;
	kin.S_x = S_x;
	kin.S_p = S_p;
	kin.V_1 = V_1;
	kin.V_2 = V_2;
	kin.V_m = V_m;
	kin.V_p = V_p;

	kin.lambda_S = lambda_S;
	kin.lambda_Y = lambda_Y;
	kin.lambda_1 = lambda_1;
	kin.lambda_2 = lambda_2;
	kin.lambda_3 = lambda_3;
	kin.lambda_S_sqrt = lambda_S_sqrt;
	kin.lambda_Y_sqrt = lambda_Y_sqrt;
	kin.lambda_1_sqrt = lambda_1_sqrt;
	kin.lambda_2_sqrt = lambda_2_sqrt;
	kin.lambda_3_sqrt = lambda_3_sqrt;

	kin.ph_0 = ph_0;
	kin.ph_t = ph_t;
	kin.ph_l = ph_l;
	kin.q_0 = q_0;
	kin.q_t = q_t;
	kin.q_l = q_l;
	kin.k2_0 = k2_0;
	kin.k2_t = k2_t;
	kin.k2_l = k2_l;
	kin.k1_t = k1_t;
	kin.mx_sq = mx_sq;
	kin.mx = mx;
	kin.vol_phi_h = vol_phi_h;

	kin.C_1 = C_1;

	return kin;
}

Kinematics KinematicsRad::project_shift() const {
	Kinematics kin;
	kin.target = target;
	kin.beam = beam;
	kin.hadron = hadron;

	kin.S = S;
	kin.M = M;
	kin.m = m;
	kin.mh = mh;
	kin.M_th = M_th;

	kin.x = shift_x;
	kin.y = shift_y;
	kin.z = shift_z;
	kin.ph_t_sq = shift_ph_t_sq;
	kin.phi_h = shift_phi_h;
	kin.phi = phi;

	kin.phi_q = shift_phi_q;

	kin.cos_phi_h = shift_cos_phi_h;
	kin.sin_phi_h = shift_sin_phi_h;
	kin.cos_phi = cos_phi;
	kin.sin_phi = sin_phi;
	kin.cos_phi_q = shift_cos_phi_q;
	kin.sin_phi_q = shift_sin_phi_q;

	kin.Q_sq = shift_Q_sq;
	kin.Q = shift_Q;
	kin.t = shift_t;
	kin.X = X;
	kin.S_x = shift_S_x;
	kin.S_p = S_p;
	kin.V_1 = V_1;
	kin.V_2 = V_2;
	kin.V_m = shift_V_m;
	kin.V_p = V_p;

	kin.lambda_S = lambda_S;
	kin.lambda_Y = shift_lambda_Y;
	kin.lambda_1 = shift_lambda_1;
	kin.lambda_2 = shift_lambda_2;
	kin.lambda_3 = shift_lambda_3;
	kin.lambda_S_sqrt = lambda_S_sqrt;
	kin.lambda_Y_sqrt = shift_lambda_Y_sqrt;
	kin.lambda_1_sqrt = shift_lambda_1_sqrt;
	kin.lambda_2_sqrt = shift_lambda_2_sqrt;
	kin.lambda_3_sqrt = shift_lambda_3_sqrt;

	kin.ph_0 = ph_0;
	kin.ph_t = shift_ph_t;
	kin.ph_l = shift_ph_l;
	kin.q_0 = shift_q_0;
	kin.q_t = shift_q_t;
	kin.q_l = shift_q_l;
	kin.k2_0 = k2_0;
	kin.k2_t = k2_t;
	kin.k2_l = k2_l;
	kin.k1_t = shift_k1_t;
	kin.mx_sq = shift_mx_sq;
	kin.mx = shift_mx;
	kin.vol_phi_h = shift_vol_phi_h;

	kin.C_1 = shift_C_1;

	return kin;
}

Final::Final(Initial init, Vec3 target_pol, Kinematics kin) {
	Transform4 target = lab_from_target(init, target_pol);
	Transform4 hadron = target * target_from_hadron(kin);
	// `q` is easy to reconstruct in the hadron frame, since the z-axis is
	// defined to point along `q`.
	q = hadron * Vec4(kin.q_0, 0., 0., kin.lambda_Y_sqrt/(2.*kin.M));
	k2 = target * Vec4(
		kin.k2_0,
		kin.k2_t * kin.sin_phi,
		kin.k2_t * kin.cos_phi,
		kin.k2_l);
	ph = hadron * Vec4(kin.ph_0, kin.ph_t, 0., kin.ph_l);
}

FinalRad::FinalRad(Initial init, Vec3 target_pol, KinematicsRad kin) {
	Transform4 target = lab_from_target(init, target_pol);
	Transform4 lepton = target * target_from_lepton(kin.project());
	q = lepton * Vec4(kin.q_0, 0., 0., kin.lambda_Y_sqrt/(2.*kin.M));
	k2 = target * Vec4(
		kin.k2_0,
		kin.k2_t * kin.sin_phi,
		kin.k2_t * kin.cos_phi,
		kin.k2_l);
	// To be slightly more efficient, construct both the `ph` and `k` vectors in
	// the lepton frame, as they are simply rotated by `phi_h` and `phi_k` about
	// the z-axis in this frame.
	ph = lepton * Vec4(
		kin.ph_0,
		kin.ph_t * kin.cos_phi_h,
		kin.ph_t * kin.sin_phi_h,
		kin.ph_l);
	k = lepton * Vec4(
		kin.k_0,
		kin.k_t * kin.cos_phi_k,
		kin.k_t * kin.sin_phi_k,
		kin.k_l);
}

// Bounds of kinematic variables.
// TODO: Some of the calculations in this section are redundant with earlier
// kinematic calculations. This should be refactored to avoid that later.
Bounds kin::x_bounds(Initial) {
	return Bounds(0., 1.);
}
Bounds kin::y_bounds(Initial init, Real x) {
	Real M = mass(init.target);
	Real m = mass(init.beam);
	Real S = 2.*dot(init.p, init.k1);
	Real lambda_S = sq(S) - 4.*sq(M)*sq(m);
	// TODO: Include bound from `sq(px) >= 0`.
	Real min = 0.;
	Real max = std::fmin(
		1.,
		// Bound from the condition that `sq(q_t) >= 0`.
		(x*lambda_S)/(S*(sq(m) + sq(M)*sq(x) + S*x)));
	return Bounds(min, max);
}
Bounds kin::z_bounds(Initial init, Hadron h, Real M_th, Real x, Real y) {
	Real M = mass(init.target);
	Real mh = mass(h);
	Real S = 2.*dot(init.p, init.k1);
	Real S_x = S*y;
	Real Q_sq = S*x*y;
	Real lambda_Y_sqrt = std::sqrt(sq(S_x) + 4.*sq(M)*Q_sq);
	Real lambda = sq(M) + S_x - Q_sq;
	Real denom = 2.*S_x*lambda;
	Real a = (S_x + 2.*sq(M))*(lambda + sq(mh) - sq(M_th));
	Real b = (lambda - sq(mh))*lambda_Y_sqrt;
	Real min = std::fmax(
		0.,
		// Bound from the condition that `mx_sq >= 0`.
		(a - b)/denom);
	Real max = std::fmin(
		1.,
		// Bound from the condition that `mx_sq >= 0`.
		(a + b)/denom);
	return Bounds(min, max);
}
Bounds kin::ph_t_sq_bounds(Initial init, Hadron h, Real M_th, Real x, Real y, Real z) {
	Real M = mass(init.target);
	Real mh = mass(h);
	Real S = 2.*dot(init.p, init.k1);
	Real S_x = S*y;
	Real Q_sq = S*x*y;
	Real lambda_Y = sq(S_x) + 4.*sq(M)*Q_sq;
	Real lambda = sq(M) + S_x - Q_sq;
	Real ph_0 = (z*S_x)/(2.*M);
	Real min = 0.;
	// Bound from the condition that `mx_sq >= 0`.
	Real max = sq(ph_0) - sq(mh) - sq(M)/lambda_Y*sq(
		lambda + sq(mh) - sq(M_th) - S_x*z - (2.*sq(ph_0))/z);
	return Bounds(min, max);
}
Bounds kin::tau_bounds(Kinematics kin) {
	// Equation [1.44].
	Real min = (kin.S_x - kin.lambda_Y_sqrt)/(2.*sq(kin.M));
	Real max = (kin.S_x + kin.lambda_Y_sqrt)/(2.*sq(kin.M));
	return Bounds(min, max);
}
Bounds kin::R_bounds(Kinematics kin, Real tau, Real phi_k) {
	Real tau_min = tau_bounds(kin).min;
	Real tau_max = tau_bounds(kin).max;
	// Copied from the kinematic calculations below.
	Real mu = kin.ph_0/kin.M
		+ 1./kin.lambda_Y_sqrt*(
			(2.*tau*sq(kin.M) - kin.S_x)*kin.ph_l/kin.M
			- 2.*kin.M*kin.ph_t*std::cos(kin.phi_h - phi_k)*std::sqrt(
				(tau - tau_min)*(tau_max - tau)));
	// Equation [1.44].
	Real min = 0.;
	Real max = 1./(1. + tau - mu)*(kin.mx_sq - sq(kin.M_th));
	return Bounds(min, max);
}
Bounds kin::R_bounds_soft(Kinematics kin, Real tau, Real phi_k, Real k0_cut) {
	Real k0_max = (kin.mx_sq - sq(kin.M_th))/(2.*kin.mx);
	Real tau_min = tau_bounds(kin).min;
	Real tau_max = tau_bounds(kin).max;
	// Copied from the kinematic calculation for `mu`.
	Real mu = kin.ph_0/kin.M
		+ 1./kin.lambda_Y_sqrt*(
			(2.*tau*sq(kin.M) - kin.S_x)*kin.ph_l/kin.M
			- 2.*kin.M*kin.ph_t*std::cos(kin.phi_h - phi_k)*std::sqrt(
				(tau - tau_min)*(tau_max - tau)));
	// Equation [1.44].
	Real min = 0.;
	Real max = 1./(1. + tau - mu)*(
		k0_cut < k0_max ?
		(2.*kin.mx)*k0_cut :
		kin.mx_sq - sq(kin.M_th));
	if (min <= max) {
		return Bounds(min, max);
	} else {
		return Bounds(min, min);
	}
}
Bounds kin::R_bounds_hard(Kinematics kin, Real tau, Real phi_k, Real k0_cut) {
	Real tau_min = tau_bounds(kin).min;
	Real tau_max = tau_bounds(kin).max;
	Real mu = kin.ph_0/kin.M
		+ 1./kin.lambda_Y_sqrt*(
			(2.*tau*sq(kin.M) - kin.S_x)*kin.ph_l/kin.M
			- 2.*kin.M*kin.ph_t*std::cos(kin.phi_h - phi_k)*std::sqrt(
				(tau - tau_min)*(tau_max - tau)));
	Real min = 1./(1. + tau - mu)*(
		k0_cut > 0. ?
		2.*kin.mx*k0_cut :
		0.);
	Real max = 1./(1. + tau - mu)*(kin.mx_sq - sq(kin.M_th));
	if (min <= max) {
		return Bounds(min, max);
	} else {
		return Bounds(max, max);
	}
}

// Check whether within kinematic bounds.
bool kin::valid(Kinematics kin) {
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
	} else if (!(kin.mx >= kin.M_th)) {
		return false;
	} else if (!(kin.mx_sq <= kin.S + sq(kin.M) - sq(kin.m) - sq(kin.mh))) {
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
bool kin::valid(KinematicsRad kin_rad) {
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

