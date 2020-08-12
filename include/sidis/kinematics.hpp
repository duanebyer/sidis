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

