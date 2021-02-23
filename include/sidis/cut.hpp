#ifndef SIDIS_CUT_HPP
#define SIDIS_CUT_HPP

#include "sidis/numeric.hpp"
#include "sidis/extra/bound.hpp"

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
	// TODO: Cuts should be possible to enable or disable. (Bound are options?)
	// Base cuts.
	math::Bound x;
	math::Bound y;
	math::Bound z;
	math::Bound ph_t_sq;
	math::Bound phi_h;
	math::Bound phi;

	// Kinematic variable cuts.
	// TODO: Special cut.
	math::Bound Q_sq;
	math::Bound t;
	math::Bound w;

	// Energy cuts.
	math::Bound mx_sq;
	math::Bound q_0;
	math::Bound k2_0;
	math::Bound ph_0;

	// Angle cuts.
	math::Bound theta_q;
	math::Bound theta_k2;
	math::Bound theta_h;

	Cut();
};

struct CutRad {
	math::Bound tau;
	math::Bound phi_k;
	math::Bound k_0_bar;

	// Energy cuts.
	math::Bound k_0;

	// Angle cuts.
	math::Bound theta_k;

	CutRad();
};

Real S_min(kin::Particles ps);

math::Bound x_bounds(kin::Particles ps, Real S);
math::Bound y_bounds(kin::Particles ps, Real S, Real x);
math::Bound z_bounds(kin::Particles ps, Real S, Real x, Real y);
math::Bound ph_t_sq_bounds(kin::Particles ps, Real S, Real x, Real y, Real z);
math::Bound tau_bounds(kin::Kinematics kin);
math::Bound R_bounds(kin::Kinematics kin, Real tau, Real phi_k);

math::Bound x_bounds(Cut cut, kin::Particles ps, Real S);
math::Bound y_bounds(Cut cut, kin::Particles ps, Real S, Real x);
math::Bound z_bounds(Cut cut, kin::Particles ps, Real S, Real x, Real y);
math::Bound ph_t_sq_bounds(Cut cut, kin::Particles ps, Real S, Real x, Real y, Real z);
math::Bound tau_bounds(CutRad cut, kin::Kinematics kin);
math::Bound R_bounds(CutRad cut, kin::Kinematics kin, Real tau, Real phi_k);

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

