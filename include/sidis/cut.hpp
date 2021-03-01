#ifndef SIDIS_CUT_HPP
#define SIDIS_CUT_HPP

#include "sidis/bound.hpp"
#include "sidis/numeric.hpp"

namespace sidis {

namespace kin {
	struct PhaseSpace;
	struct PhaseSpaceRad;
	struct Kinematics;
	struct KinematicsRad;
}

namespace part {
	struct Particles;
}

namespace cut {

struct Cut {
	// Base cuts.
	math::Bound x;
	math::Bound y;
	math::Bound z;
	math::Bound ph_t_sq;
	math::Bound phi_h;
	math::Bound phi;

	// Kinematic variable cuts.
	math::Bound Q_sq;
	math::Bound t;
	math::Bound w;
	math::Bound r;

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

Real S_min(part::Particles const& ps);

math::Bound x_bound(part::Particles const& ps, Real S);
math::Bound y_bound(part::Particles const& ps, Real S, Real x);
math::Bound z_bound(part::Particles const& ps, Real S, Real x, Real y);
math::Bound ph_t_sq_bound(part::Particles const& ps, Real S, Real x, Real y, Real z);
math::Bound tau_bound(kin::Kinematics const& kin);
math::Bound R_bound(kin::Kinematics const& kin, Real tau, Real phi_k);

math::Bound x_bound(Cut const& cut, part::Particles const& ps, Real S);
math::Bound y_bound(Cut const& cut, part::Particles const& ps, Real S, Real x);
math::Bound z_bound(Cut const& cut, part::Particles const& ps, Real S, Real x, Real y);
math::Bound ph_t_sq_bound(Cut const& cut, part::Particles const& ps, Real S, Real x, Real y, Real z);
math::Bound tau_bound(CutRad const& cut, kin::Kinematics const& kin);
math::Bound R_bound(CutRad const& cut, kin::Kinematics const& kin, Real tau, Real phi_k);

bool valid(kin::Kinematics const& kin);
bool valid(kin::KinematicsRad const& kin);
bool valid(Cut const& cut, kin::Kinematics const& kin);
bool valid(Cut const& cut, CutRad const& cut_rad, kin::KinematicsRad const& kin);
bool valid(CutRad const& cut, kin::KinematicsRad const& kin);

// Draw from the non-radiative kinematic region.
bool take(
	part::Particles const& ps, Real S, const Real point[6],
	kin::PhaseSpace* ph_space_out, Real* jacobian_out);
bool take(
	part::Particles const& ps, Real S, const Real point[6],
	kin::Kinematics* kin_out, Real* jacobian_out);
bool take(
	Cut const& cut,
	part::Particles const& ps, Real S, const Real point[6],
	kin::Kinematics* kin_out, Real* jacobian_out);

// Draw from the radiative kinematic region.
bool take(
	part::Particles const& ps, Real S, const Real point[9],
	kin::PhaseSpaceRad* ph_space_out, Real* jacobian_out);
bool take(
	part::Particles const& ps, Real S, const Real point[9],
	kin::KinematicsRad* kin_out, Real* jacobian_out);
bool take(
	Cut const& cut, CutRad const& cut_rad,
	part::Particles const& ps, Real S, const Real point[9],
	kin::KinematicsRad* kin_out, Real* jacobian_out);
bool take(
	kin::Kinematics const& kin, const Real point[3],
	kin::PhaseSpaceRad* ph_space_out, Real* jacobian_out);
bool take(
	kin::Kinematics const& kin, const Real point[3],
	kin::KinematicsRad* kin_out, Real* jacobian_out);
bool take(
	CutRad const& cut, kin::Kinematics const& kin, const Real point[3],
	kin::KinematicsRad* kin_out, Real* jacobian_out);

}
}

#endif

