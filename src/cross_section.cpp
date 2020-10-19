#include "sidis/cross_section.hpp"

#include <cmath>
#include <limits>
#include <stdexcept>

#include <iostream>
#include <ios>
#include <cubature.hpp>

#include "sidis/constant.hpp"
#include "sidis/frame.hpp"
#include "sidis/hadronic_coeff.hpp"
#include "sidis/kinematics.hpp"
#include "sidis/leptonic_coeff.hpp"
#include "sidis/structure_function.hpp"
#include "sidis/extra/math.hpp"
#include "sidis/extra/integrate.hpp"
#include "sidis/extra/transform.hpp"
#include "sidis/extra/vector.hpp"

using namespace sidis;
using namespace sidis::xs;
using namespace sidis::constant;
using namespace sidis::had;
using namespace sidis::integ;
using namespace sidis::kin;
using namespace sidis::lep;
using namespace sidis::math;
using namespace sidis::sf;

namespace {

Real const SMALL_R_REL = std::cbrt(2. * std::numeric_limits<Real>::epsilon());
Real const DELTA_R_REL = std::sqrt(2. * std::numeric_limits<Real>::epsilon());

Real S_phi(Real z_il[4], Real z_ij[4][4], Real a, Real b) {
	// Equation [1.40].
	return -a*(
		b*std::log((z_il[0]*z_il[2])/(z_il[1]*z_il[3]))
		+ 0.5*(
			sq(std::log(std::abs(z_il[0])))
			- sq(std::log(std::abs(z_il[1])))
			- sq(std::log(std::abs(z_il[2])))
			+ sq(std::log(std::abs(z_il[3]))))
		+ (
			std::log(std::abs(z_il[0]))*std::log(std::abs(z_ij[0][1]))
			- std::log(std::abs(z_il[0]))*std::log(std::abs(z_ij[0][2]))
			- std::log(std::abs(z_il[0]))*std::log(std::abs(z_ij[0][3]))
			- std::log(std::abs(z_il[1]))*std::log(std::abs(z_ij[1][0]))
			+ std::log(std::abs(z_il[1]))*std::log(std::abs(z_ij[1][2]))
			+ std::log(std::abs(z_il[1]))*std::log(std::abs(z_ij[1][3]))
			+ std::log(std::abs(z_il[2]))*std::log(std::abs(z_ij[2][0]))
			+ std::log(std::abs(z_il[2]))*std::log(std::abs(z_ij[2][1]))
			- std::log(std::abs(z_il[2]))*std::log(std::abs(z_ij[2][3]))
			- std::log(std::abs(z_il[3]))*std::log(std::abs(z_ij[3][0]))
			- std::log(std::abs(z_il[3]))*std::log(std::abs(z_ij[3][1]))
			+ std::log(std::abs(z_il[3]))*std::log(std::abs(z_ij[3][2])))
		- (
			dilog(z_il[0]/z_ij[0][1])
			- dilog(z_il[0]/z_ij[0][2])
			- dilog(z_il[0]/z_ij[0][3])
			- dilog(z_il[1]/z_ij[1][0])
			+ dilog(z_il[1]/z_ij[1][2])
			+ dilog(z_il[1]/z_ij[1][3])
			+ dilog(z_il[2]/z_ij[2][0])
			+ dilog(z_il[2]/z_ij[2][1])
			- dilog(z_il[2]/z_ij[2][3])
			- dilog(z_il[3]/z_ij[3][0])
			- dilog(z_il[3]/z_ij[3][1])
			+ dilog(z_il[3]/z_ij[3][2])));
}

}

Born::Born(Kinematics kin) {
	// Equation [1.15]. The `Q^4` factor has been absorbed into `C_1`.
	coeff = (sq(ALPHA)*kin.S*sq(kin.S_x))
		/(8.*kin.M*kin.ph_l*kin.lambda_S);
}

Amm::Amm(Kinematics kin) {
	// Equation [1.53]. The `Q^4` factor has been absorbed into `C_1`.
	Real lambda_m = kin.Q_sq*(kin.Q_sq + 4.*sq(kin.m));
	Real lambda_m_sqrt = std::sqrt(lambda_m);
	Real diff_m = sqrt1p_1m((4.*sq(kin.m))/kin.Q_sq);
	Real sum_m = 2. + diff_m;
	Real L_m = 1./lambda_m_sqrt*std::log(sum_m/diff_m);
	coeff = L_m*kin.Q_sq*(std::pow(ALPHA, 3)*sq(kin.m)*kin.S*sq(kin.S_x))
		/(16.*PI*kin.M*kin.ph_l*kin.lambda_S);
}

Rad::Rad(KinematicsRad kin) {
	// Equation [1.43].
	coeff = -(std::pow(ALPHA, 3)*kin.S*sq(kin.S_x))
		/(64.*sq(PI)*kin.M*kin.ph_l*kin.lambda_S*kin.lambda_Y_sqrt);
	R = kin.R;
}

Real xs::born(Real lambda_e, Vec3 eta, Kinematics kin, Model const& model) {
	if (model.target != kin.target) {
		throw std::invalid_argument(
			"Structure functions are for wrong target nucleus");
	}
	Born b(kin);
	SfXX sf = model.sf(kin.hadron, kin.x, kin.z, kin.Q_sq, kin.ph_t_sq);
	LepBornXX lep(kin);
	HadXX had(kin, sf);

	Real uu = born_uu(b, lep, had);
	Vec3 up(
		born_ut1(b, lep, had),
		born_ut2(b, lep, had),
		born_ul(b, lep, had));
	Real lu = born_lu(b, lep, had);
	Vec3 lp(
		born_lt1(b, lep, had),
		born_lt2(b, lep, had),
		born_ll(b, lep, had));
	return uu + dot(eta, up) + lambda_e * (lu + dot(eta, lp));
}

Real xs::amm(Real lambda_e, Vec3 eta, Kinematics kin, Model const& model) {
	if (model.target != kin.target) {
		throw std::invalid_argument(
			"Structure functions are for wrong target nucleus");
	}

	Amm b(kin);
	SfXX sf = model.sf(kin.hadron, kin.x, kin.z, kin.Q_sq, kin.ph_t_sq);
	LepAmmXX lep(kin);
	HadXX had(kin, sf);

	Real uu = amm_uu(b, lep, had);
	Vec3 up(
		amm_ut1(b, lep, had),
		amm_ut2(b, lep, had),
		amm_ul(b, lep, had));
	Real lu = amm_lu(b, lep, had);
	Vec3 lp(
		amm_lt1(b, lep, had),
		amm_lt2(b, lep, had),
		amm_ll(b, lep, had));
	return uu + dot(eta, up) + lambda_e * (lu + dot(eta, lp));
}

Real xs::nrad(Real lambda_e, Vec3 eta, Kinematics kin, Model const& model) {
	return xs::born_rad_factor(kin) * xs::born(lambda_e, eta, kin, model)
		+ xs::amm(lambda_e, eta, kin, model);
}

Real xs::rad(Real lambda_e, Vec3 eta, KinematicsRad kin, Model const& model) {
	return rad_hard(lambda_e, eta, kin, model)
		+ rad_soft(lambda_e, eta, kin, model);
}

Real xs::rad_integ(Real lambda_e, math::Vec3 eta, kin::Kinematics kin, sf::Model const& model) {
	cubature::EstErr<Real> xs_integ = cubature::cubature<3>(
		[&](cubature::Point<3, Real> x) {
			Bounds tau_b = tau_bounds(kin);
			Real tau = tau_b.lerp(x[0]);

			Bounds phi_k_b(0., 2. * PI);
			Real phi_k = phi_k_b.lerp(x[1]);

			Bounds R_b = R_bounds(kin, tau, phi_k);
			Real R = R_b.lerp(x[2]);

			Real jacobian = tau_b.size() * phi_k_b.size() * R_b.size();
			KinematicsRad kin_rad(kin, tau, phi_k, R);
			Real result = rad(lambda_e, eta, kin_rad, model);
			if (std::isnan(result)) {
				// If the result is `NaN`, it most likely means we went out of
				// the allowed region for the structure function grids. In that
				// case, just return zero.
				// TODO: Handle this case in a more correct way.
				return 0.;
			} else {
				return jacobian * result;
			}
		},
		cubature::Point<3, Real>{ 0., 0., 0. },
		cubature::Point<3, Real>{ 1., 1., 1. },
		10000, 0., 1e-6);
	return xs_integ.val;
}

Real xs::rad_hard(Real lambda_e, Vec3 eta, KinematicsRad kin, Model const& model) {
	if (model.target != kin.target) {
		throw std::invalid_argument(
			"Structure functions are for wrong target nucleus");
	}

	Rad b(kin);
	// A shifted version of `eta` is needed to properly combine the unshifted
	// and shifted cross-sections.
	Transform3 shift_rot = frame::hadron_from_shift(kin);
	Sf shift_sf = model.sf(
		kin.hadron,
		kin.shift_x, kin.shift_z, kin.shift_Q_sq, kin.shift_ph_t_sq);

	LepRadXX lep(kin);
	HadXX shift_had(kin.project_shift(), shift_sf);

	Real uu = rad_hard_uu(b, lep, shift_had);
	Vec3 up = rad_hard_up(b, lep, shift_had, shift_rot);
	Real lu = rad_hard_lu(b, lep, shift_had);
	Vec3 lp = rad_hard_lp(b, lep, shift_had, shift_rot);
	return uu + dot(eta, up) + lambda_e * (lu + dot(eta, lp));
}

Real xs::rad_soft(Real lambda_e, Vec3 eta, KinematicsRad kin, Model const& model) {
	if (model.target != kin.target) {
		throw std::invalid_argument(
			"Structure functions are for wrong target nucleus");
	}
	if (std::abs(kin.R) < std::abs(kin.R_max) * SMALL_R_REL) {
		// For small `R`, the evaluation of the soft contribution to the
		// radiative cross-section becomes inaccurate due to catastrophic
		// cancellation.
		return rad_soft_0(lambda_e, eta, kin, model);
	} else {
		Rad b(kin);
		Transform3 shift_rot = frame::hadron_from_shift(kin);
		Sf sf = model.sf(
			kin.hadron,
			kin.x, kin.z, kin.Q_sq, kin.ph_t_sq);
		Sf shift_sf = model.sf(
			kin.hadron,
			kin.shift_x, kin.shift_z, kin.shift_Q_sq, kin.shift_ph_t_sq);
		LepRadXX lep(kin);
		HadXX had(kin.project(), sf);
		HadXX shift_had(kin.project_shift(), shift_sf);

		Real uu = rad_soft_uu_R(b, lep, shift_had);
		Vec3 up = rad_soft_up_R(b, lep, shift_had, shift_rot);
		Real lu = rad_soft_lu_R(b, lep, shift_had);
		Vec3 lp = rad_soft_lp_R(b, lep, shift_had, shift_rot);

		Real uu_0 = rad_soft_uu_R0(b, lep, had);
		Vec3 up_0(
			rad_soft_ul_R0(b, lep, had),
			rad_soft_ut1_R0(b, lep, had),
			rad_soft_ut2_R0(b, lep, had));
		Real lu_0 = rad_soft_lu_R0(b, lep, had);
		Vec3 lp_0(
			rad_soft_ll_R0(b, lep, had),
			rad_soft_lt1_R0(b, lep, had),
			rad_soft_lt2_R0(b, lep, had));
		Real xs_R = uu + dot(eta, up) + lambda_e * (lu + dot(eta, lp));
		Real xs_R0 = uu_0 + dot(eta, up_0) + lambda_e * (lu_0 + dot(eta, lp_0));
		return 1. / kin.R * (xs_R - xs_R0);
	}
}

Real xs::rad_soft_0(Real lambda_e, Vec3 eta, KinematicsRad kin, Model const& model) {
	if (model.target != kin.target) {
		throw std::invalid_argument(
			"Structure functions are for wrong target nucleus");
	}
	// For small `R`, the evaluation of the soft contribution to the radiative
	// cross-section becomes inaccurate due to catastrophic cancellation. At the
	// same time, this is the region that contributes the most to the cross-
	// section, so it's important to get right. We use linear extrapolation to
	// extend the radiative cross-section to this region. The first and second
	// derivatives of the cross-section are used for this.
	Real delta_R = kin.R_max * DELTA_R_REL;
	Real R_rel = kin.R / delta_R;
	Real R_1 = delta_R;
	Real R_2 = 2. * delta_R;
	Rad b(kin);
	KinematicsRad kin_1(kin.project(), kin.tau, kin.phi_k, R_1);
	KinematicsRad kin_2(kin.project(), kin.tau, kin.phi_k, R_2);
	Transform3 shift_rot_1 = frame::hadron_from_shift(kin_1);
	Transform3 shift_rot_2 = frame::hadron_from_shift(kin_2);
	Sf sf_0 = model.sf(
		kin_1.hadron,
		kin_1.x, kin_1.z, kin_1.Q_sq, kin_1.ph_t_sq);
	Sf sf_1 = model.sf(
		kin_1.hadron,
		kin_1.shift_x, kin_1.shift_z, kin_1.shift_Q_sq, kin_1.shift_ph_t_sq);
	Sf sf_2 = model.sf(
		kin_2.hadron,
		kin_2.shift_x, kin_2.shift_z, kin_2.shift_Q_sq, kin_2.shift_ph_t_sq);
	LepRadXX lep(kin);
	HadXX had_0(kin.project(), sf_0);
	HadXX had_1(kin_1.project_shift(), sf_1);
	HadXX had_2(kin_2.project_shift(), sf_2);
	Real uu_0 = rad_soft_uu_R0(b, lep, had_0);
	Vec3 up_0(
		rad_soft_ul_R0(b, lep, had_0),
		rad_soft_ut1_R0(b, lep, had_0),
		rad_soft_ut2_R0(b, lep, had_0));
	Real lu_0 = rad_soft_lu_R0(b, lep, had_0);
	Vec3 lp_0(
		rad_soft_ll_R0(b, lep, had_0),
		rad_soft_lt1_R0(b, lep, had_0),
		rad_soft_lt2_R0(b, lep, had_0));
	Real uu_1 = rad_soft_uu_R(b, lep, had_1);
	Vec3 up_1 = rad_soft_up_R(b, lep, had_1, shift_rot_1);
	Real lu_1 = rad_soft_lu_R(b, lep, had_1);
	Vec3 lp_1 = rad_soft_lp_R(b, lep, had_1, shift_rot_1);
	Real uu_2 = rad_soft_uu_R(b, lep, had_2);
	Vec3 up_2 = rad_soft_up_R(b, lep, had_2, shift_rot_2);
	Real lu_2 = rad_soft_lu_R(b, lep, had_2);
	Vec3 lp_2 = rad_soft_lp_R(b, lep, had_2, shift_rot_2);
	Real xs_R0 = uu_0 + dot(eta, up_0) + lambda_e * (lu_0 + dot(eta, lp_0));
	Real xs_R1 = uu_1 + dot(eta, up_1) + lambda_e * (lu_1 + dot(eta, lp_1));
	Real xs_R2 = uu_2 + dot(eta, up_2) + lambda_e * (lu_2 + dot(eta, lp_2));

	// Compute first and second derivatives, and use the result as part of a
	// linear extrapolation for the cross-section.
	return 1. / delta_R * (
		(xs_R1 - xs_R0)
		+ 0.5 * R_rel * (xs_R2 - 2. * xs_R1 + xs_R0));
}

// Born base functions.
Real xs::born_uu(Born b, LepBornUU lep, HadUU had) {
	return b.coeff*(
		lep.theta_1*had.H_10
		+ lep.theta_2*had.H_20
		+ lep.theta_3*had.H_30
		+ lep.theta_4*had.H_40);
}
Real xs::born_ul(Born b, LepBornUP lep, HadUL had) {
	return b.coeff*(lep.theta_6*had.H_63 + lep.theta_8*had.H_83);
}
Real xs::born_ut1(Born b, LepBornUP lep, HadUT1 had) {
	return b.coeff*(lep.theta_6*had.H_61 + lep.theta_8*had.H_81);
}
Real xs::born_ut2(Born b, LepBornUU lep, HadUT2 had) {
	return b.coeff*(
		lep.theta_1*had.H_12
		+ lep.theta_2*had.H_22
		+ lep.theta_3*had.H_32
		+ lep.theta_4*had.H_42);
}
Real xs::born_lu(Born b, LepBornLU lep, HadLU had) {
	return b.coeff*lep.theta_5*had.H_50;
}
Real xs::born_ll(Born b, LepBornLP lep, HadLL had) {
	return b.coeff*(lep.theta_7*had.H_73 + lep.theta_9*had.H_93);
}
Real xs::born_lt1(Born b, LepBornLP lep, HadLT1 had) {
	return b.coeff*(lep.theta_7*had.H_71 + lep.theta_9*had.H_91);
}
Real xs::born_lt2(Born b, LepBornLU lep, HadLT2 had) {
	return b.coeff*lep.theta_5*had.H_52;
}

// AMM base functions.
Real xs::amm_uu(Amm b, LepAmmUU lep, HadUU had) {
	return b.coeff*(
		lep.theta_1*had.H_10
		+ lep.theta_2*had.H_20
		+ lep.theta_3*had.H_30
		+ lep.theta_4*had.H_40);
}
Real xs::amm_ul(Amm b, LepAmmUP lep, HadUL had) {
	return b.coeff*(lep.theta_6*had.H_63 + lep.theta_8*had.H_83);
}
Real xs::amm_ut1(Amm b, LepAmmUP lep, HadUT1 had) {
	return b.coeff*(lep.theta_6*had.H_61 + lep.theta_8*had.H_81);
}
Real xs::amm_ut2(Amm b, LepAmmUU lep, HadUT2 had) {
	return b.coeff*(
		lep.theta_1*had.H_12
		+ lep.theta_2*had.H_22
		+ lep.theta_3*had.H_32
		+ lep.theta_4*had.H_42);
}
Real xs::amm_lu(Amm b, LepAmmLU lep, HadLU had) {
	return b.coeff*lep.theta_5*had.H_50;
}
Real xs::amm_ll(Amm b, LepAmmLP lep, HadLL had) {
	return b.coeff*(lep.theta_7*had.H_73 + lep.theta_9*had.H_93);
}
Real xs::amm_lt1(Amm b, LepAmmLP lep, HadLT1 had) {
	return b.coeff*(lep.theta_7*had.H_71 + lep.theta_9*had.H_91);
}
Real xs::amm_lt2(Amm b, LepAmmLU lep, HadLT2 had) {
	return b.coeff*lep.theta_5*had.H_52;
}

// Hard radiative base functions.
Real xs::rad_hard_uu(Rad b, LepRadUU lep, HadUU shift_had) {
	return b.coeff*(
		(
			lep.theta_012*shift_had.H_10
			+ lep.theta_022*shift_had.H_20
			+ lep.theta_032*shift_had.H_30
			+ lep.theta_042*shift_had.H_40)
		+ b.R*(
			lep.theta_013*shift_had.H_10
			+ lep.theta_023*shift_had.H_20
			+ lep.theta_033*shift_had.H_30
			+ lep.theta_043*shift_had.H_40));
}
Vec3 xs::rad_hard_up(Rad b, LepRadUX lep, HadUP shift_had, Transform3 shift_rot) {
	return b.coeff*dot(shift_rot, Vec3(
		// UT1.
		(lep.theta_062*shift_had.H_61 + lep.theta_082*shift_had.H_81)
		+ b.R*(lep.theta_063*shift_had.H_61 + lep.theta_083*shift_had.H_81)
		+ b.R*b.R*(lep.theta_064*shift_had.H_61 + lep.theta_084*shift_had.H_81),
		// UT2.
		(
			lep.theta_012*shift_had.H_12
			+ lep.theta_022*shift_had.H_22
			+ lep.theta_032*shift_had.H_32
			+ lep.theta_042*shift_had.H_42)
		+ b.R*(
			lep.theta_013*shift_had.H_12
			+ lep.theta_023*shift_had.H_22
			+ lep.theta_033*shift_had.H_32
			+ lep.theta_043*shift_had.H_42),
		// UL.
		(lep.theta_062*shift_had.H_63 + lep.theta_082*shift_had.H_83)
		+ b.R*(lep.theta_063*shift_had.H_63 + lep.theta_083*shift_had.H_83)
		+ b.R*b.R*(lep.theta_064*shift_had.H_63 + lep.theta_084*shift_had.H_83)));
}
Real xs::rad_hard_lu(Rad b, LepRadLU lep, HadLU shift_had) {
	return b.coeff*(
		(lep.theta_052 + lep.theta_152)*shift_had.H_50
		+ b.R*(lep.theta_053 + lep.theta_153)*shift_had.H_50);
}
Vec3 xs::rad_hard_lp(Rad b, LepRadLX lep, HadLP shift_had, Transform3 shift_rot) {
	return b.coeff*dot(shift_rot, Vec3(
		// LT1.
		(
			(lep.theta_072 + lep.theta_172)*shift_had.H_71
			+ (lep.theta_092 + lep.theta_192)*shift_had.H_91)
		+ b.R*(
			(lep.theta_073 + lep.theta_173)*shift_had.H_71
			+ (lep.theta_093 + lep.theta_193)*shift_had.H_91)
		+ b.R*b.R*(
			(lep.theta_074 + lep.theta_174)*shift_had.H_71
			+ (lep.theta_094 + lep.theta_194)*shift_had.H_91),
		// LT2.
		(lep.theta_052 + lep.theta_152)*shift_had.H_52
		+ b.R*(lep.theta_053 + lep.theta_153)*shift_had.H_52,
		// LL.
		(
			(lep.theta_072 + lep.theta_172)*shift_had.H_73
			+ (lep.theta_092 + lep.theta_192)*shift_had.H_93)
		+ b.R*(
			(lep.theta_073 + lep.theta_173)*shift_had.H_73
			+ (lep.theta_093 + lep.theta_193)*shift_had.H_93)
		+ b.R*b.R*(
			(lep.theta_074 + lep.theta_174)*shift_had.H_73
			+ (lep.theta_094 + lep.theta_194)*shift_had.H_93)));
}

// Soft radiative base functions.
Real xs::rad_soft_uu_R(Rad b, LepRadUU lep, HadUU shift_had) {
	return b.coeff*(
			lep.theta_011*shift_had.H_10
			+ lep.theta_021*shift_had.H_20
			+ lep.theta_031*shift_had.H_30
			+ lep.theta_041*shift_had.H_40);
}
Vec3 xs::rad_soft_up_R(Rad b, LepRadUX lep, HadUP shift_had, Transform3 shift_rot) {
	return b.coeff*dot(shift_rot, Vec3(
		// UT1.
		lep.theta_061*shift_had.H_61 + lep.theta_081*shift_had.H_81,
		// UT2.
		lep.theta_011*shift_had.H_12
		+ lep.theta_021*shift_had.H_22
		+ lep.theta_031*shift_had.H_32
		+ lep.theta_041*shift_had.H_42,
		// UL.
		lep.theta_061*shift_had.H_63 + lep.theta_081*shift_had.H_83));
}
Real xs::rad_soft_lu_R(Rad b, LepRadLU lep, HadLU shift_had) {
	return b.coeff*(lep.theta_051 + lep.theta_151)*shift_had.H_50;
}
Vec3 xs::rad_soft_lp_R(Rad b, LepRadLX lep, HadLP shift_had, Transform3 shift_rot) {
	return b.coeff*dot(shift_rot, Vec3(
		// LT1.
		(lep.theta_071 + lep.theta_171)*shift_had.H_71
		+ (lep.theta_091 + lep.theta_191)*shift_had.H_91,
		// LT2.
		(lep.theta_051 + lep.theta_151)*shift_had.H_52,
		// LL.
		(lep.theta_071 + lep.theta_171)*shift_had.H_73
		+ (lep.theta_091 + lep.theta_191)*shift_had.H_93));
}

// Soft radiative base functions with hadronic part at `R=0`.
Real xs::rad_soft_uu_R0(Rad b, LepRadUU lep, HadUU had) {
	return rad_soft_uu_R(b, lep, had);
}
Real xs::rad_soft_ul_R0(Rad b, LepRadUP lep, HadUL had) {
	return b.coeff*(lep.theta_061*had.H_63 + lep.theta_081*had.H_83);
}
Real xs::rad_soft_ut1_R0(Rad b, LepRadUP lep, HadUT1 had) {
	return b.coeff*(lep.theta_061*had.H_61 + lep.theta_081*had.H_81);
}
Real xs::rad_soft_ut2_R0(Rad b, LepRadUU lep, HadUT2 had) {
	return b.coeff*(
		lep.theta_011*had.H_12
		+ lep.theta_021*had.H_22
		+ lep.theta_031*had.H_32
		+ lep.theta_041*had.H_42);
}
Real xs::rad_soft_lu_R0(Rad b, LepRadLU lep, HadLU had) {
	return rad_soft_lu_R(b, lep, had);
}
Real xs::rad_soft_ll_R0(Rad b, LepRadLP lep, HadLL had) {
	return b.coeff*(
		(lep.theta_071 + lep.theta_171)*had.H_73
		+ (lep.theta_091 + lep.theta_191)*had.H_93);
}
Real xs::rad_soft_lt1_R0(Rad b, LepRadLP lep, HadLT1 had) {
	return b.coeff*(
		(lep.theta_071 + lep.theta_171)*had.H_71
		+ (lep.theta_091 + lep.theta_191)*had.H_91);
}
Real xs::rad_soft_lt2_R0(Rad b, LepRadLU lep, HadLT2 had) {
	return b.coeff*(lep.theta_051 + lep.theta_151)*had.H_52;
}

// Radiative correction to Born cross-section.
Real xs::born_rad_factor(Kinematics kin) {
	Real vr = delta_vr(kin);
	Real had = delta_vac_had(kin);
	Real lep = delta_vac_lep(kin);
	return 1. + ALPHA/PI*(vr + had + lep);
}

Real xs::delta_vr(Kinematics kin) {
	// Equation [1.3].
	Real Q_m_sq = kin.Q_sq + 2.*sq(kin.m);
	Real S_prime = kin.S - kin.Q_sq - kin.V_1;
	Real X_prime = kin.X + kin.Q_sq - kin.V_2;
	Real lambda_m = kin.Q_sq*(kin.Q_sq + 4.*sq(kin.m));
	Real lambda_S_prime = sq(S_prime) - 4.*sq(kin.m)*kin.mx_sq;
	Real lambda_X_prime = sq(X_prime) - 4.*sq(kin.m)*kin.mx_sq;
	Real lambda_m_sqrt = std::sqrt(lambda_m);
	Real lambda_S_prime_sqrt = std::sqrt(lambda_S_prime);
	Real lambda_X_prime_sqrt = std::sqrt(lambda_X_prime);

	// Differences of the form `√λ/|S| - 1`.
	Real diff_m = sqrt1p_1m((4.*sq(kin.m))/kin.Q_sq);
	Real diff_S_prime = sqrt1p_1m(-(4.*sq(kin.m)*kin.mx_sq)/sq(S_prime));
	Real diff_X_prime = sqrt1p_1m(-(4.*sq(kin.m)*kin.mx_sq)/sq(X_prime));
	Real sum_m = 2. + diff_m;
	Real sum_S_prime = 2. + diff_S_prime;
	Real sum_X_prime = 2. + diff_X_prime;

	// Equation [1.C10].
	Real L_m = 1./lambda_m_sqrt*std::log(sum_m/diff_m);
	Real L_S_prime = 1./lambda_S_prime_sqrt*std::log(-sum_S_prime/diff_S_prime);
	Real L_X_prime = 1./lambda_X_prime_sqrt*std::log(-sum_X_prime/diff_X_prime);

	// Instead of computing `z_i` using equation [1.42], compute the differences
	// directly. Use the notation `z_ij = z_j - z_i`.
	Real c = (2.*kin.mx_sq*kin.Q_sq)/X_prime;
	Real c_diff = c/diff_X_prime;
	Real c_sum = c/sum_X_prime;
	Real d = lambda_X_prime_sqrt*diff_m;
	Real d_diff = X_prime*diff_X_prime;
	Real d_sum = X_prime*sum_X_prime;

	Real z_1u = 1./lambda_X_prime_sqrt*(
		S_prime*sum_S_prime - X_prime*sum_X_prime - c_diff*diff_m);
	Real z_2u = 1./lambda_X_prime_sqrt*(
		S_prime*sum_S_prime - X_prime*sum_X_prime + c_diff*sum_m);
	Real z_3u = 1./lambda_X_prime_sqrt*(
		S_prime*diff_S_prime - X_prime*diff_X_prime + c_sum*sum_m);
	Real z_4u = 1./lambda_X_prime_sqrt*(
		S_prime*diff_S_prime - X_prime*diff_X_prime - c_sum*diff_m);
	Real z_1d = 1./lambda_X_prime*(
		X_prime*sum_X_prime*(S_prime - X_prime) - c_diff*(d + d_diff));
	Real z_2d = 1./lambda_X_prime*(
		X_prime*sum_X_prime*(S_prime - X_prime) + c_diff*(d + d_sum));
	Real z_3d = 1./lambda_X_prime*(
		-X_prime*diff_X_prime*(S_prime - X_prime) + c_sum*(d + d_diff));
	Real z_4d = 1./lambda_X_prime*(
		-X_prime*diff_X_prime*(S_prime - X_prime) - c_sum*(d + d_sum));

	Real z_12 = -lambda_m_sqrt/lambda_X_prime_sqrt
		*(4.*kin.mx_sq)/(X_prime*diff_X_prime);
	Real z_34 = lambda_m_sqrt/lambda_X_prime_sqrt
		*(4.*kin.mx_sq)/(X_prime*sum_X_prime);
	Real z_13 = 1./lambda_X_prime_sqrt*(
		2.*(S_prime - X_prime) - c*(sum_m/sum_X_prime + diff_m/diff_X_prime));
	Real z_24 = 1./lambda_X_prime_sqrt*(
		2.*(S_prime - X_prime) + c*(diff_m/sum_X_prime + sum_m/diff_X_prime));
	Real z_14 = 1./lambda_X_prime_sqrt*(
		2.*(S_prime - X_prime) - 2.*c*diff_m/(diff_X_prime*sum_X_prime));
	Real z_23 = 1./lambda_X_prime_sqrt*(
		2.*(S_prime - X_prime) + 2.*c*sum_m/(diff_X_prime*sum_X_prime));

	Real z_iu[4] = { z_1u, z_2u, z_3u, z_4u };
	Real z_id[4] = { z_1d, z_2d, z_3d, z_4d };
	Real z_ij[4][4] = {
		{    0.,  z_12,  z_13,  z_14 },
		{ -z_12,    0.,  z_23,  z_24 },
		{ -z_13, -z_23,    0.,  z_34 },
		{ -z_14, -z_24, -z_34,    0. },
	};

	// Equation [1.40].
	Real S_phi_a = Q_m_sq/(2.*lambda_m_sqrt);
	Real S_phi_b = std::log(-diff_X_prime/sum_X_prime);
	Real S_phi_diff = S_phi(z_iu, z_ij, S_phi_a, S_phi_b)
		- S_phi(z_id, z_ij, S_phi_a, S_phi_b);

	// Equation [1.52].
	Real delta =
		2.*(Q_m_sq*L_m - 1.)
			*std::log((kin.mx_sq - sq(kin.M_th))/(kin.m*kin.mx))
		+ 0.5*S_prime*L_S_prime + 0.5*X_prime*L_X_prime
		+ S_phi_diff
		- 2.
		+ (1.5*kin.Q_sq + 4.*sq(kin.m))*L_m
		- Q_m_sq/lambda_m_sqrt*(
			0.5*lambda_m*sq(L_m)
			+ 2.*dilog((2.*lambda_m_sqrt)/(kin.Q_sq + lambda_m_sqrt))
			- 0.5*sq(PI));
	return delta;
}

Real xs::delta_vac_lep(Kinematics kin) {
	// Equation [1.50].
	Real ms[3] = { MASS_E, MASS_MU, MASS_TAU };
	Real delta = 0.;
	for (unsigned idx = 0; idx < 3; ++idx) {
		Real m = ms[idx];
		Real lambda_sqrt = std::sqrt(kin.Q_sq*(kin.Q_sq + 4.*sq(m)));
		Real diff_m = sqrt1p_1m((4.*sq(m))/kin.Q_sq);
		Real sum_m = 2. + diff_m;
		Real L_m = 1/lambda_sqrt*std::log(sum_m/diff_m);
		delta += 2./3.*(kin.Q_sq
			+ 2.*sq(m))*L_m
			- 10./9.
			+ (8.*sq(m))/(3.*kin.Q_sq)*(1. - 2.*sq(m)*L_m);
	}
	return delta;
}

Real xs::delta_vac_had(Kinematics kin) {
	if (kin.Q_sq < 1.) {
		return -(2.*PI)/ALPHA*(-1.345e-9 - 2.302e-3*std::log(1. + 4.091*kin.Q_sq));
	} else if (kin.Q_sq < 64.) {
		return -(2.*PI)/ALPHA*(-1.512e-3 - 2.822e-3*std::log(1. + 1.218*kin.Q_sq));
	} else {
		return -(2.*PI)/ALPHA*(-1.1344-3 - 3.0680-3*std::log(1. + 0.99992*kin.Q_sq));
	}
}

