#ifndef SIDIS_CUT_HPP
#define SIDIS_CUT_HPP

#include "sidis/numeric.hpp"
#include "sidis/extra/bounds.hpp"

namespace sidis {

namespace kin {
	struct Particles;
	struct PhaseSpace;
	struct PhaseSpaceRad;
	struct Kinematics;
	struct KinematicsRad;
}

namespace cut {

struct Cut {
	math::Bounds x;
	math::Bounds y;
	math::Bounds z;
	math::Bounds ph_t_sq;
	math::Bounds phi_h;
	math::Bounds phi;

	Cut();
};

struct CutRad {
	math::Bounds tau;
	math::Bounds phi_k;
	math::Bounds k_0_bar;

	CutRad();
};

Real S_min(kin::Particles ps);

math::Bounds x_bounds(kin::Particles ps, Real S);
math::Bounds y_bounds(kin::Particles ps, Real S, Real x);
math::Bounds z_bounds(kin::Particles ps, Real S, Real x, Real y);
math::Bounds ph_t_sq_bounds(kin::Particles ps, Real S, Real x, Real y, Real z);
math::Bounds tau_bounds(kin::Kinematics kin);
math::Bounds R_bounds(kin::Kinematics kin, Real tau, Real phi_k);

math::Bounds x_bounds(Cut cut, kin::Particles ps, Real S);
math::Bounds y_bounds(Cut cut, kin::Particles ps, Real S, Real x);
math::Bounds z_bounds(Cut cut, kin::Particles ps, Real S, Real x, Real y);
math::Bounds ph_t_sq_bounds(Cut cut, kin::Particles ps, Real S, Real x, Real y, Real z);
math::Bounds tau_bounds(CutRad cut, kin::Kinematics kin);
math::Bounds R_bounds(CutRad cut, kin::Kinematics kin, Real tau, Real phi_k);

bool valid(kin::Kinematics kin);
bool valid(kin::KinematicsRad kin);
bool valid(Cut cut, kin::Kinematics kin);
bool valid(Cut cut, CutRad cut_rad, kin::KinematicsRad kin);
bool valid(CutRad cut, kin::KinematicsRad kin);

// Draw from the non-radiative kinematic region.
bool take(
	kin::Particles ps, Real S, const Real point[6],
	kin::PhaseSpace* ph_space_out, Real* jacobian_out);
bool take(
	kin::Particles ps, Real S, const Real point[6],
	kin::Kinematics* kin_out, Real* jacobian_out);
bool take(
	Cut cut, kin::Particles ps, Real S, const Real point[6],
	kin::Kinematics* kin_out, Real* jacobian_out);

// Draw from the radiative kinematic region.
bool take(
	kin::Particles ps, Real S, const Real point[9],
	kin::PhaseSpaceRad* ph_space_out, Real* jacobian_out);
bool take(
	kin::Particles ps, Real S, const Real point[9],
	kin::KinematicsRad* kin_out, Real* jacobian_out);
bool take(
	Cut cut, CutRad cut_rad, kin::Particles ps, Real S, const Real point[9],
	kin::KinematicsRad* kin_out, Real* jacobian_out);
bool take(
	kin::Kinematics kin, const Real point[3],
	kin::PhaseSpaceRad* ph_space_out, Real* jacobian_out);
bool take(
	kin::Kinematics kin, const Real point[3],
	kin::KinematicsRad* kin_out, Real* jacobian_out);
bool take(
	CutRad cut, kin::Kinematics kin, const Real point[3],
	kin::KinematicsRad* kin_out, Real* jacobian_out);

}
}

#endif

