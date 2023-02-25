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

/**
 * \defgroup CutGroup Cuts
 *
 * Definitions of kinematic cuts on the SIDIS process.
 */
/// \{

/**
 * Contains all supported cuts that can be applied to the non-radiative SIDIS
 * process. If any of the cuts are invalid, then that cut is ignored.
 */
struct Cut {
	/// Base cuts
	/// \{
	math::Bound x;
	math::Bound y;
	math::Bound z;
	math::Bound ph_t_sq;
	math::Bound phi_h;
	math::Bound phi;
	/// \}

	/// Additional cuts
	/// \{
	math::Bound Q_sq;
	math::Bound t;
	math::Bound W_sq;
	math::Bound r;
	math::Bound mx_sq;
	math::Bound qt_to_Q;
	/// \}

	/// Momentum cuts in lab frame
	/// \{
	math::Bound lab_mom_q;
	math::Bound lab_mom_k2;
	math::Bound lab_mom_h;
	/// \}

	/// Angle cuts in lab frame
	/// \{
	math::Bound lab_theta_q;
	math::Bound lab_theta_k2;
	math::Bound lab_theta_h;
	/// \}
};

/**
 * Contains additional cuts for use with radiative SIDIS processes only.
 */
struct CutRad {
	/// Base cuts
	/// \{
	math::Bound tau;
	math::Bound phi_k;
	math::Bound R;
	math::Bound k_0_bar;
	/// \}

	/// Momentum cuts in lab frame
	/// \{
	math::Bound lab_mom_k;
	math::Bound lab_theta_k;
	/// \}
};

/// Minimum value that \f$S\f$ can take for a given set of particles.
Real S_min(part::Particles const& ps);

/// \name Kinematic limits
/// These functions give a set of kinematic limits on the six base SIDIS
/// variables. These limits are less restrictive than the true SIDIS kinematic
/// limits, so all points outside of these limits are guaranteed invalid, but a
/// point inside may or may not be valid. To check validity precisely, see the
/// \ref cut::valid() functions.
///
/// The six base SIDIS variables are chosen one at a time, in the order \f$x\f$,
/// \f$y\f$, \f$z\f$, \f$p_t^2\f$. The angles can be chosen at any time, and can
/// have any value. For each variable that is fixed to a certain value, new
/// bounds are then imposed on the remaining variables that are yet to be
/// chosen. The following methods provide an approximation of these further
/// bounds.
///
/// For radiative SIDIS processes, there are additional restrictions on
/// \f$\tau\f$ and \f$R\f$, in that order.
/// \{

/// Kinematic limits on \f$x\f$.
math::Bound x_bound(part::Particles const& ps, Real S);
/// Kinematic limits on \f$y\f$.
math::Bound y_bound(part::Particles const& ps, Real S, Real x);
/// Kinematic limits on \f$z\f$.
math::Bound z_bound(part::Particles const& ps, Real S, Real x, Real y);
/// Kinematic limits on \f$p_t^2\f$.
math::Bound ph_t_sq_bound(part::Particles const& ps, Real S, Real x, Real y, Real z);
/// Kinematic limits on \f$\tau\f$.
math::Bound tau_bound(kin::Kinematics const& kin);
/// Kinematic limits on \f$R\f$.
math::Bound R_bound(kin::Kinematics const& kin, Real tau, Real phi_k);
/// \}

/// \name Cut limits
/// These methods are similar to the kinematic limits (see methods like
/// x_bound()), but additionally take into account any user-defined cuts on
/// various kinematic variables. Remember that it is necessary to use the
/// cut::valid() methods to verify that a point chosen within these bounds is
/// really within the desired phase space region.
/// \{

/// Kinematic limits on \f$x\f$, including Cut%s.
math::Bound x_bound(Cut const& cut, part::Particles const& ps, Real S);
/// Kinematic limits on \f$y\f$, including Cut%s.
math::Bound y_bound(Cut const& cut, part::Particles const& ps, Real S, Real x);
/// Kinematic limits on \f$z\f$, including Cut%s.
math::Bound z_bound(Cut const& cut, part::Particles const& ps, Real S, Real x, Real y);
/// Kinematic limits on \f$p_t^2\f$, including Cut%s.
math::Bound ph_t_sq_bound(Cut const& cut, part::Particles const& ps, Real S, Real x, Real y, Real z);
/// Kinematic limits on \f$\tau\f$, including Cut%s.
math::Bound tau_bound(CutRad const& cut, kin::Kinematics const& kin);
/// Kinematic limits on \f$R\f$, including Cut%s.
math::Bound R_bound(CutRad const& cut, kin::Kinematics const& kin, Real tau, Real phi_k);
/// \}

/// \name Validity checks
/// These functions check that the kinematics of a SIDIS process are
/// kinematically valid, and optionally, also whether they satisfy
/// user-provided cuts.
/// \{

/// Check whether a set of Kinematics are within the allowed phase space of
/// non-radiative SIDIS.
bool valid(kin::Kinematics const& kin);
/// Check whether a set of KinematicsRad are within the allowed phase space of
/// radiative SIDIS.
bool valid(kin::KinematicsRad const& kin);
/// Check whether a set of Kinematics satisfy some kinematic Cut%s. Also checks
/// kinematic validity.
bool valid(Cut const& cut, kin::Kinematics const& kin);
/// Check whether a set of KinematicsRad satisfy some kinematic Cut%s. Also
/// checks kinematic validity.
bool valid(Cut const& cut, CutRad const& cut_rad, kin::KinematicsRad const& kin);
/// Check whether a set of KinematicsRad satisfy some kinematic Cut%s. Also
/// checks kinematic validity.
bool valid(CutRad const& cut, kin::KinematicsRad const& kin);
/// \}

/// \name Phase space sampling
/// These functions draw from the valid phase space using a point in the unit
/// hyper-cube. Returns whether they were successful in sampling or not. It's
/// recommended to use these methods in combination with a random number
/// generator.
/// \{

/// Draw from non-radiative SIDIS phase space.
bool take(
	part::Particles const& ps, Real S, const Real point[6],
	kin::PhaseSpace* ph_space_out, Real* jac_out);
/// Draw from non-radiative SIDIS phase space.
bool take(
	part::Particles const& ps, Real S, const Real point[6],
	kin::Kinematics* kin_out, Real* jac_out);
/// Draw from non-radiative SIDIS phase space subject to Cut%s.
bool take(
	Cut const& cut,
	part::Particles const& ps, Real S, const Real point[6],
	kin::Kinematics* kin_out, Real* jac_out);

/// Draw from the radiative SIDIS phase space.
bool take(
	part::Particles const& ps, Real S, const Real point[9],
	kin::PhaseSpaceRad* ph_space_out, Real* jac_out);
/// Draw from the radiative SIDIS phase space.
bool take(
	part::Particles const& ps, Real S, const Real point[9],
	kin::KinematicsRad* kin_out, Real* jac_out);
/// Draw from the radiative SIDIS phase space subject to Cut%s.
bool take(
	Cut const& cut, CutRad const& cut_rad,
	part::Particles const& ps, Real S, const Real point[9],
	kin::KinematicsRad* kin_out, Real* jac_out);
/// Draw from the radiative SIDIS phase space.
bool take(
	kin::Kinematics const& kin, const Real point[3],
	kin::PhaseSpaceRad* ph_space_out, Real* jac_out);
/// Draw from the radiative SIDIS phase space.
bool take(
	kin::Kinematics const& kin, const Real point[3],
	kin::KinematicsRad* kin_out, Real* jac_out);
/// Draw from the radiative SIDIS phase space subject to Cut%s.
bool take(
	CutRad const& cut, kin::Kinematics const& kin, const Real point[3],
	kin::KinematicsRad* kin_out, Real* jac_out);
/// \}
/// \}

}
}

#endif

