#ifndef __SIDIS_KINEMATICS_HPP__
#define __SIDIS_KINEMATICS_HPP__

#include "sidis/numeric.hpp"
#include "sidis/vector.hpp"

namespace sidis {
namespace kin {

struct Initial {
	math::Vec4 p;
	math::Vec4 k1;
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
	Real k_t;
	Real mx_sq;
	Real mx;
	Real vol_phi_h;

	Kinematics(Initial init, PhaseSpace ph_space, Real mh, Real M_th);
};

struct KinematicsRad {
};

struct KinematicsEx {
};

struct Final {
	math::Vec4 q;
	math::Vec4 k2;
	math::Vec4 ph;
	math::Vec4 px;

	Final(Initial init, Kinematics kin);
};

struct FinalRad {
	math::Vec4 q;
	math::Vec4 k2;
	math::Vec4 ph;
	math::Vec4 px;
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

