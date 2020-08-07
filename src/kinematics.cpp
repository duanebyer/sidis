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
	M = init.p.norm();
	m = init.k1.norm();

	// Equation [1.3].
	Q_sq = S*x*y;
	Q = std::sqrt(Q_sq);
	X = S*(1. - y);
	S_x = S - X;
	S_p = S + X;
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
	ph_l = ph_0*std::sqrt(1. - ph_t_sq/sq(ph_0) - sq(mh/ph_0));
	t = (2.*M*ph_l*lambda_Y_sqrt - z*sq(S_x)) / (2.*sq(M)) - Q_sq + sq(mh);
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

Final::Final(Initial init, Kinematics kin) {
	// First reconstruct in the target rest frame, then transform into the same
	// frame as was used for the initial state.
	Real phi_h = kin.phi_h;
	Real phi = kin.phi;
	Real M = kin.M;
	Real Q_sq = kin.Q_sq;
	Real S = kin.S;
	Real S_x = kin.S_x;
	Real ph_0 = kin.ph_0;
	Real ph_l = kin.ph_l;
	Real ph_t = kin.ph_t;

	// Equations [A.1], [A.2].
	Real q_0 = S_x/(2.*M);
	Real q_l = (2.*sq(M)*Q_sq + S*S_x)/(2.*M*kin.lambda_S_sqrt);
	Real q_t = q_0*std::sqrt(1. - sq(q_l/q_0) + Q_sq/sq(q_0));
	q = Vec4(q_0, -q_t*std::cos(phi), -q_t*std::sin(phi), q_l);

	k2 = init.k1 - q;

	// Form a basis for the reconstruction of ph.
	Vec3 eq_y = cross(q.r, init.k1.r).unit();
	Vec3 eq_z = q.r.unit();
	Vec3 eq_x = cross(eq_y, eq_z);

	ph = Vec4(
		ph_0,
		ph_l*eq_z + ph_t*std::cos(phi_h)*eq_x + ph_t*std::sin(phi_h)*eq_y);
}

