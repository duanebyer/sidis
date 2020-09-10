#ifndef SIDIS_KINEMATICS_HPP
#define SIDIS_KINEMATICS_HPP

#include "sidis/numeric.hpp"
#include "sidis/vector.hpp"

namespace sidis {
namespace kin {

struct Initial {
	Real M;
	Real m;
	math::Vec4 p;
	math::Vec4 k1;

	Initial(Real M, math::Vec3 p, Real m, math::Vec3 k1) :
		M(M),
		m(m),
		p(math::Vec4::from_length_and_r(M, p)),
		k1(math::Vec4::from_length_and_r(m, k1)) { }

	Initial(Real M, Real m, Real beam_energy) :
		M(M),
		m(m),
		p(math::Vec4(M, 0., 0., 0.)),
		k1(math::Vec4::from_length_and_t(m, beam_energy, math::Vec3::Z)) { }
};

struct PhaseSpace {
	Real x;
	Real y;
	Real z;
	Real ph_t_sq;
	Real phi_h;
	Real phi;
};

struct PhaseSpaceRad {
	Real x;
	Real y;
	Real z;
	Real ph_t_sq;
	Real phi_h;
	Real phi;

	Real tau;
	Real phi_k;
	Real R;
};

struct PhaseSpaceEx {
	Real x;
	Real y;
	Real z;
	Real ph_t_sq;
	Real phi_h;
	Real phi;

	Real tau;
	Real phi_k;
};

struct Kinematics {
	Real S;
	Real M;
	Real m;
	Real mh;
	Real M_th;

	Real x;
	Real y;
	Real z;
	Real ph_t_sq;
	Real phi_h;
	Real phi;

	Real Q_sq;
	Real Q;
	Real t;
	Real X;
	Real S_x;
	Real S_p;
	Real V_1;
	Real V_2;
	Real V_m;
	Real V_p;

	Real lambda_S;
	Real lambda_Y;
	Real lambda_1;
	Real lambda_2;
	Real lambda_3;
	Real lambda_S_sqrt;
	Real lambda_Y_sqrt;
	Real lambda_1_sqrt;
	Real lambda_2_sqrt;
	Real lambda_3_sqrt;

	Real ph_0;
	Real ph_t;
	Real ph_l;
	Real k1_t;
	Real mx_sq;
	Real mx;
	Real vol_phi_h;

	Kinematics() { }
	Kinematics(Initial init, PhaseSpace ph_space, Real mh, Real M_th);
};

struct KinematicsRad {
	// Base set of kinematic variables.
	Real S;
	Real M;
	Real m;
	Real mh;
	Real M_th;

	Real x;
	Real y;
	Real z;
	Real ph_t_sq;
	Real phi_h;
	Real phi;

	Real Q_sq;
	Real Q;
	Real t;
	Real X;
	Real S_x;
	Real S_p;
	Real V_1;
	Real V_2;
	Real V_m;
	Real V_p;

	Real lambda_S;
	Real lambda_Y;
	Real lambda_1;
	Real lambda_2;
	Real lambda_3;
	Real lambda_V;
	Real lambda_RV;
	Real lambda_RY;
	Real lambda_S_sqrt;
	Real lambda_Y_sqrt;
	Real lambda_1_sqrt;
	Real lambda_2_sqrt;
	Real lambda_3_sqrt;

	Real ph_0;
	Real ph_t;
	Real ph_l;
	Real k1_t;
	Real mx_sq;
	Real mx;
	Real vol_phi_h;

	// Additional radiative kinematic variables.
	Real tau;
	Real phi_k;
	Real R;

	Real tau_min;
	Real tau_max;
	Real R_max;
	Real mu;
	Real z_1;
	Real z_2;

	Real lambda_z;
	Real lambda_H;
	Real lambda_z_sqrt;

	Real vol_phi_k;
	Real vol_phi_hk;

	Real F_22;
	Real F_21;
	Real F_2p;
	Real F_2m;
	Real F_d;
	Real F_1p;
	Real F_IR;

	// Shifted kinematic variables.
	Real shift_x;
	Real shift_y;
	Real shift_z;
	Real shift_ph_t_sq;
	Real shift_phi_h;

	Real shift_Q_sq;
	Real shift_Q;
	Real shift_t;
	Real shift_S_x;
	Real shift_V_m;

	Real shift_lambda_Y;
	Real shift_lambda_1;
	Real shift_lambda_2;
	Real shift_lambda_3;
	Real shift_lambda_Y_sqrt;
	Real shift_lambda_1_sqrt;
	Real shift_lambda_2_sqrt;
	Real shift_lambda_3_sqrt;

	Real shift_ph_t;
	Real shift_ph_l;
	Real shift_k1_t;
	Real shift_mx_sq;
	Real shift_mx;
	Real shift_vol_phi_h;

	Kinematics project() const;
	Kinematics project_shift() const;

	KinematicsRad() { }
	KinematicsRad(
		Initial init,
		PhaseSpaceRad ph_space,
		Real mh, Real M_th) :
		KinematicsRad(
			Kinematics(
				init,
				{
					ph_space.x,
					ph_space.y,
					ph_space.z,
					ph_space.ph_t_sq,
					ph_space.phi_h,
					ph_space.phi,
				},
				mh, M_th),
			ph_space.tau,
			ph_space.phi_k,
			ph_space.R) { }
	KinematicsRad(Kinematics kin, Real tau, Real phi_k, Real R);
};

struct KinematicsEx {
};

struct Final {
	math::Vec4 q;
	math::Vec4 k2;
	math::Vec4 ph;

	Final(Initial init, Kinematics kin);
};

struct FinalRad {
	math::Vec4 q;
	math::Vec4 k2;
	math::Vec4 ph;
	math::Vec4 k;

	FinalRad(Initial init, KinematicsRad kin);
};

struct FinalEx {
	math::Vec4 q;
	math::Vec4 k2;
	math::Vec4 ph;
	math::Vec4 pu;
	math::Vec4 k;

	FinalEx(Initial init, KinematicsEx kin);
};

}
}

#endif

