#include "sidis/kinematics.hpp"

#include <cmath>

#include "sidis/math.hpp"

using namespace sidis;
using namespace sidis::kin;
using namespace sidis::math;

Kinematics::Kinematics(Initial init, PhaseSpace ph_space, Real mh, Real M_th) :
		mh(mh),
		M_th(M_th) {
	x = ph_space.x;
	y = ph_space.y;
	z = ph_space.z;
	ph_t_sq = ph_space.ph_t_sq;
	phi_h = ph_space.phi_h;
	phi = ph_space.phi;

	S = 2. * dot(init.p, init.k1);
	M = init.M;
	m = init.m;

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
	// In the low `ph_t` case (where the cross-section is the highest), the
	// computation for `t` has a catastrophic cancellation between the terms
	// `2 M ph_l √λ_Y - z S_x²`. So, it's better to compute `t` in the following
	// way:
	Real lambda_Y_ratio = (4.*sq(M)*Q_sq)/sq(S_x);
	t = -Q_sq + sq(mh) + (ph_0*S_x)/M*sqrt1p_1m(
		lambda_Y_ratio - ph_ratio_sq - lambda_Y_ratio*ph_ratio_sq);
	mx_sq = sq(M) + t + (1. - z)*S_x;
	mx = std::sqrt(mx_sq);
	k_t = lambda_1_sqrt/lambda_Y_sqrt;

	// Equation [1.5].
	V_1 = ph_0*S/M - (ph_l*(S*S_x + 2.*sq(M)*Q_sq))/(M*lambda_Y_sqrt)
		- 2.*ph_t*k_t*std::cos(phi_h);
	V_2 = ph_0*X/M - (ph_l*(X*S_x - 2.*sq(M)*Q_sq))/(M*lambda_Y_sqrt)
		- 2.*ph_t*k_t*std::cos(phi_h);
	V_p = 0.5*(V_1 + V_2);
	V_m = 0.5*(sq(mh) - Q_sq - t);

	// Paragraph below equation [1.14].
	lambda_2 = sq(V_m) + sq(mh)*Q_sq;
	lambda_3 = V_m + z*Q_sq;
	lambda_2_sqrt = std::sqrt(lambda_2);
	lambda_3_sqrt = std::sqrt(lambda_3);

	// Equation [1.6]. `vol_phi_h` is defined as `dot(epsilon_perp, ph)`.
	vol_phi_h = -0.5*ph_t*lambda_1_sqrt*std::sin(phi_h);
}

KinematicsRad::KinematicsRad(Kinematics kin, Real tau, Real phi_k, Real R) :
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
		k_t(kin.k_t),
		mx_sq(kin.mx_sq),
		mx(kin.mx),
		vol_phi_h(kin.vol_phi_h),
		tau(tau),
		phi_k(phi_k),
		R(R) {
	// TODO: Fill in equation number from derivations.
	lambda_H = sq(z*S_x) - 4.*sq(M)*sq(mh);
	// TODO: Fill in equation number from derivations.
	lambda_V = z*sq(S_x) - 4.*sq(M)*V_m;
	lambda_RY = R*(S_x - 2.*sq(M)*tau);
	lambda_RV = (2.*M)/lambda_Y_sqrt*(
		std::sqrt((sq(R)*lambda_Y - sq(lambda_RY))*ph_t_sq)
			*std::cos(phi_h - phi_k)
		+ lambda_RY*ph_l);

	// Equation [1.44].
	tau_min = (S_x - lambda_Y_sqrt)/(2.*sq(M));
	tau_max = (S_x + lambda_Y_sqrt)/(2.*sq(M));

	// Equation [1.B3].
	mu = ph_0/M
		+ (ph_l*(2.*tau*sq(M) - S_x))/(M*lambda_Y_sqrt)
		- 2.*M*ph_t*std::cos(phi_h + phi_k)
			*std::sqrt((tau_max - tau)*(tau - tau_min))/lambda_Y_sqrt;

	// Equation [1.44].
	R_max = (mx_sq - sq(M_th))/(1. + tau - mu);

	// Equation [1.B4].
	lambda_z = (tau_max - tau)*(tau - tau_min)*lambda_1;
	lambda_z_sqrt = std::sqrt(lambda_z);
	z_1 = 1./lambda_Y*(
		Q_sq*S_p
		+ tau*(S*S_x + 2.*sq(M)*Q_sq)
		- 2.*M*lambda_z_sqrt*std::cos(phi_k));
	z_2 = 1./lambda_Y*(
		Q_sq*S_p
		+ tau*(X*S_x - 2.*sq(M)*Q_sq)
		- 2.*M*lambda_z_sqrt*std::cos(phi_k));

	// Equation [1.B5].
	F_22 = 1./sq(z_2);
	F_21 = 1./sq(z_1);
	F_2p = F_22 + F_21;
	F_2m = F_22 - F_21;
	F_d = 1./(z_1*z_2);
	F_1p = 1./z_1 + 1./z_2;
	F_IR = sq(m)*F_2p - (Q_sq + 2.*sq(m))*F_d;

	// Equation [1.30].
	// TODO: Why does this equation require a negative sign compared to what is
	// given in [1]?
	vol_phi_k = -(std::sin(phi_k)*R*std::sqrt(
			lambda_1*(Q_sq + tau*(S_x - tau*sq(M)))))
		/(2.*lambda_Y_sqrt);
	// Equation [1.A9].
	vol_phi_hk = 1./(2.*lambda_1)*(
		R*vol_phi_h*(z_1*lambda_Y - Q_sq*S_p - tau*(S*S_x + 2.*sq(M)*Q_sq))
		+ vol_phi_k*(S_x*(z*Q_sq*S_p - S*V_2 + X*V_1) - 4.*V_p*sq(M)*Q_sq));

	shift_Q_sq = Q_sq + R*tau;
	shift_Q = std::sqrt(shift_Q_sq);
	shift_S_x = S_x - R;
	shift_V_m = V_m - lambda_RV/(4.*sq(M));

	shift_x = shift_Q_sq/shift_S_x;
	shift_z = (2.*M*ph_0)/shift_S_x;
	// TODO: Fill in equation number in derivations.
	shift_t = t - R*tau + (z*R*S_x - lambda_RV)/(2.*sq(M));

	shift_lambda_Y = lambda_Y + sq(R) - 2.*lambda_RY;
	shift_lambda_1 = shift_Q_sq*(S*X - sq(M)*shift_Q_sq) - sq(m)*shift_lambda_Y;
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
	shift_k_t = shift_lambda_1_sqrt / shift_lambda_Y_sqrt;
	shift_mx_sq = mx_sq - R*(1. + tau) + (z*R*S_x - lambda_RV)/(2.*sq(M));
	shift_mx = std::sqrt(shift_mx_sq);

	// Equation [1.A9].
	shift_vol_phi_h = vol_phi_h - 1./(2.*lambda_1)*(
		R*vol_phi_h*(
			tau*lambda_S
			+ 2.*sq(m)*S_x
			+ Q_sq*S
			- z_1*(S*S_x + 2.*sq(M)*Q_sq))
		+ vol_phi_k*(
			2.*sq(m)*(4.*V_m*sq(M) - z*sq(S_x))
			+ S*(S*V_2 - X*V_1 - z*Q_sq*S_x)
			+ 2.*V_1*sq(M)*Q_sq));

	// TODO: Fill in equation number from derivations.
	shift_phi_h = std::atan2(
		8.*sq(M)*shift_lambda_Y_sqrt*shift_vol_phi_h,
		S*shift_S_x*(shift_z*shift_lambda_Y - lambda_V)
			+ 2.*sq(M)*(lambda_V*shift_Q_sq - shift_lambda_Y*V_1));
}

Kinematics KinematicsRad::project() const {
	Kinematics kin;
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
	kin.k_t = k_t;
	kin.mx_sq = mx_sq;
	kin.mx = mx;
	kin.vol_phi_h = vol_phi_h;

	return kin;
}

Kinematics KinematicsRad::project_shift() const {
	Kinematics kin;
	kin.S = S;
	kin.M = M;
	kin.m = m;
	kin.mh = mh;
	kin.M_th = M_th;

	kin.x = shift_x;
	kin.y = y;
	kin.z = shift_z;
	kin.ph_t_sq = shift_ph_t_sq;
	kin.phi_h = shift_phi_h;
	kin.phi = phi;

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
	kin.k_t = shift_k_t;
	kin.mx_sq = shift_mx_sq;
	kin.mx = shift_mx;
	kin.vol_phi_h = shift_vol_phi_h;

	return kin;
}

Final::Final(Initial init, Kinematics kin) {
	// First reconstruct in the target rest frame, then transform into the same
	// frame as was used for the initial state.

	// Equations [A.1], [A.2].
	Real q_0 = kin.S_x/(2.*kin.M);
	Real q_l = (2.*sq(kin.M)*kin.Q_sq + kin.S*kin.S_x)
		/(2.*kin.M*kin.lambda_S_sqrt);
	Real q_t = q_0*std::sqrt(1. - sq(q_l/q_0) + kin.Q_sq/sq(q_0));
	q = Vec4(q_0, -q_t*std::cos(kin.phi), -q_t*std::sin(kin.phi), q_l);

	k2 = init.k1 - q;

	// Form a basis for the reconstruction of `ph`.
	Vec3 eq_y = cross(q.r, init.k1.r).unit();
	Vec3 eq_z = q.r.unit();
	Vec3 eq_x = cross(eq_y, eq_z);

	ph = Vec4(
		kin.ph_0,
		kin.ph_l*eq_z + kin.ph_t*(std::cos(kin.phi_h)*eq_x + std::sin(kin.phi_h)*eq_y));
}

FinalRad::FinalRad(Initial init, KinematicsRad kin) {
	Real q_0 = kin.S_x/(2.*kin.M);
	Real q_l = (2.*sq(kin.M)*kin.Q_sq + kin.S*kin.S_x)
		/(2.*kin.M*kin.lambda_S_sqrt);
	Real q_t = q_0*std::sqrt(1. - sq(q_l/q_0) + kin.Q_sq/sq(q_0));
	q = Vec4(q_0, -q_t*std::cos(kin.phi), -q_t*std::sin(kin.phi), q_l);

	k2 = init.k1 - q;

	// Both `ph` and `k` can be reconstructed with the same basis.
	Vec3 eq_y = cross(q.r, init.k1.r).unit();
	Vec3 eq_z = q.r.unit();
	Vec3 eq_x = cross(eq_y, eq_z);

	ph = Vec4(
		kin.ph_0,
		kin.ph_l*eq_z
			+ kin.ph_t*(std::cos(kin.phi_h)*eq_x + std::sin(kin.phi_h)*eq_y));

	Real k_0 = kin.R/(2.*kin.M);
	Real k_l = kin.lambda_RY/(2.*kin.M*kin.lambda_Y_sqrt);
	Real k_t = 1./(2.*kin.M)
		*std::sqrt((sq(kin.R)*kin.lambda_Y - sq(kin.lambda_RY))/kin.lambda_Y);
	k = Vec4(
		k_0,
		k_l*eq_z + k_t*(std::cos(kin.phi_k)*eq_x + std::sin(kin.phi_k)*eq_y));
}

