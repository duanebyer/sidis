#ifndef SIDIS_KINEMATICS_HPP
#define SIDIS_KINEMATICS_HPP

#include "sidis/constant.hpp"
#include "sidis/numeric.hpp"
#include "sidis/extra/vector.hpp"

namespace sidis {

namespace math {
	struct Bounds;
}

namespace kin {

struct Initial {
	constant::Nucleus target;
	constant::Lepton beam;
	math::Vec4 p;
	math::Vec4 k1;

	Initial(
		constant::Nucleus target, math::Vec3 p,
		constant::Lepton beam, math::Vec3 k1) :
		target(target),
		beam(beam),
		p(math::Vec4::from_length_and_r(mass(target), p)),
		k1(math::Vec4::from_length_and_r(mass(beam), k1)) { }

	Initial(
		constant::Nucleus target,
		constant::Lepton beam,
		Real beam_energy) :
		target(target),
		beam(beam),
		p(math::Vec4(mass(target), 0., 0., 0.)),
		k1(math::Vec4::from_length_and_t(mass(beam), beam_energy, math::Vec3::Z)) { }
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

	PhaseSpace project() const {
		return { x, y, z, ph_t_sq, phi_h, phi };
	}
};

struct Kinematics {
	constant::Nucleus target;
	constant::Lepton beam;
	constant::Hadron hadron;

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

	Real phi_q;

	Real cos_phi_h;
	Real sin_phi_h;
	Real cos_phi;
	Real sin_phi;
	Real cos_phi_q;
	Real sin_phi_q;

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

	// Components of hadron 4-momentum in lepton frame.
	Real ph_0;
	Real ph_t;
	Real ph_l;
	// Components of virtual photon 4-momentum in target frame.
	Real q_0;
	Real q_t;
	Real q_l;
	// Components of scattered lepton in target frame.
	Real k2_0;
	Real k2_t;
	Real k2_l;

	Real k1_t;
	Real mx_sq;
	Real mx;
	Real vol_phi_h;

	Real C_1;

	Kinematics() { }
	Kinematics(
		Initial init,
		PhaseSpace ph_space,
		constant::Hadron hadron,
		Real M_th);
};

struct KinematicsRad {
	// Base set of kinematic variables.
	constant::Nucleus target;
	constant::Lepton beam;
	constant::Hadron hadron;

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

	Real phi_q;

	Real cos_phi_h;
	Real sin_phi_h;
	Real cos_phi;
	Real sin_phi;
	Real cos_phi_q;
	Real sin_phi_q;

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
	Real q_0;
	Real q_t;
	Real q_l;
	Real k2_0;
	Real k2_t;
	Real k2_l;

	Real k1_t;
	Real mx_sq;
	Real mx;
	Real vol_phi_h;

	Real C_1;

	// Additional radiative kinematic variables.
	Real tau;
	Real phi_k;
	Real R;

	Real cos_phi_k;
	Real sin_phi_k;

	Real tau_min;
	Real tau_max;
	Real R_max;
	Real mu;
	Real z_1;
	Real z_2;

	Real lambda_z;
	Real lambda_H;
	Real lambda_z_sqrt;

	// Components of real photon 4-momentum in lepton frame.
	Real k_0;
	Real k_t;
	Real k_l;

	Real vol_phi_k_R;
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
	Real shift_phi_q;

	Real shift_cos_phi_h;
	Real shift_sin_phi_h;
	Real shift_cos_phi_q;
	Real shift_sin_phi_q;

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
	Real shift_q_0;
	Real shift_q_t;
	Real shift_q_l;
	Real shift_k1_t;
	Real shift_mx_sq;
	Real shift_mx;
	Real shift_vol_phi_h;

	Real shift_C_1;

	Kinematics project() const;
	Kinematics project_shift() const;

	KinematicsRad() { }
	KinematicsRad(
		Initial init,
		PhaseSpaceRad ph_space,
		constant::Hadron hadron,
		Real M_th) :
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
				hadron,
				M_th),
			ph_space.tau,
			ph_space.phi_k,
			ph_space.R) { }
	KinematicsRad(Kinematics kin, Real tau, Real phi_k, Real R);
};

struct Final {
	math::Vec4 q;
	math::Vec4 k2;
	math::Vec4 ph;

	Final(Initial init, math::Vec3 target_pol, Kinematics kin);
};

struct FinalRad {
	math::Vec4 q;
	math::Vec4 k2;
	math::Vec4 ph;
	math::Vec4 k;

	FinalRad(Initial init, math::Vec3 target_pol, KinematicsRad kin);
};

// TODO: Make an easier method for calculating the bounds.
math::Bounds x_bounds(Initial init);
math::Bounds y_bounds(Initial init, Real x);
math::Bounds z_bounds(Initial init, constant::Hadron h, Real M_th, Real x, Real y);
math::Bounds ph_t_sq_bounds(Initial init, constant::Hadron h, Real M_th, Real x, Real y, Real z);
math::Bounds tau_bounds(Kinematics kin);
math::Bounds R_bounds(Kinematics kin, Real tau, Real phi_k);
math::Bounds R_bounds_soft(Kinematics kin, Real tau, Real phi_k, Real k0_cut = constant::INF);
math::Bounds R_bounds_hard(Kinematics kin, Real tau, Real phi_k, Real k0_cut = constant::INF);

}
}

#endif

