#include "sidis/hadronic_coeff.hpp"

#include <cmath>
#include <limits>
#include <stdexcept>
#include <string>

#include "sidis/frame.hpp"
#include "sidis/kinematics.hpp"
#include "sidis/structure_function.hpp"
#include "sidis/transform.hpp"
#include "sidis/extra/exception.hpp"
#include "sidis/extra/math.hpp"

using namespace sidis;
using namespace sidis::had;
using namespace sidis::kin;
using namespace sidis::math;
using namespace sidis::sf;

namespace {

Real const DELTA_R_REL = std::sqrt(2. * std::numeric_limits<Real>::epsilon());
Real diff_R(Real R, Real delta_R, Real f_0, Real f_1, Real f_2) {
	Real R_rel = R/delta_R;
	return 1./(2.*delta_R)*((R_rel - 3.)*f_0 - 2.*(R_rel - 2.)*f_1 + (R_rel - 1.)*f_2);
}
Vec3 diff_R(Real R, Real delta_R, Vec3 f_0, Vec3 f_1, Vec3 f_2) {
	Real R_rel = R/delta_R;
	return 1./(2.*delta_R)*((R_rel - 3.)*f_0 - 2.*(R_rel - 2.)*f_1 + (R_rel - 1.)*f_2);
}

}

// Standard hadronic coefficients. Equation [1.14].
HadUU::HadUU(Kinematics kin, SfUU sf) {
	Real H_00 = kin.C_1*sf.F_UUL;
	Real H_01 = -kin.C_1*sf.F_UU_cos_phih;
	Real H_11 = kin.C_1*(sf.F_UU_cos_2phih + sf.F_UUT);
	Real H_22 = kin.C_1*(sf.F_UUT - sf.F_UU_cos_2phih);
	H_10 = H_22;
	H_20 = 4./(sq(kin.lambda_Y)*kin.ph_t_sq)*(
		kin.lambda_Y*kin.ph_t_sq*kin.Q_sq*H_00
		+ sq(kin.lambda_3)*sq(kin.S_x)*H_11
		- kin.lambda_2*kin.lambda_Y*H_22
		-2.*kin.S_x*kin.lambda_3*kin.ph_t*kin.Q*kin.lambda_Y_sqrt*H_01);
	H_30 = 1./kin.ph_t_sq*(H_11 - H_22);
	H_40 = 2./(kin.lambda_Y*kin.ph_t_sq)*(
		kin.lambda_3*kin.S_x*(H_22 - H_11)
		+ kin.ph_t*kin.Q*kin.lambda_Y_sqrt*H_01);
}
HadUL::HadUL(Kinematics kin, SfUL sf) {
	Real H_023 = kin.C_1*sf.F_UL_sin_phih;
	Real H_123 = -kin.C_1*sf.F_UL_sin_2phih;
	H_63 = 4./(std::pow(kin.lambda_Y_sqrt, 3)*kin.ph_t_sq)*(
		kin.Q*kin.ph_t*kin.lambda_Y_sqrt*H_023
		- kin.lambda_3*kin.S_x*H_123);
	H_83 = 2./(kin.lambda_Y_sqrt*kin.ph_t_sq)*H_123;
}
HadUT::HadUT(Kinematics kin, SfUT sf) {
	Real H_002 = kin.C_1*sf.F_UTL_sin_phih_m_phis;
	Real H_012 = kin.C_1*(sf.F_UT_sin_phis - sf.F_UT_sin_2phih_m_phis);
	Real H_021 = kin.C_1*(sf.F_UT_sin_2phih_m_phis + sf.F_UT_sin_phis);
	Real H_112 = kin.C_1*(sf.F_UT_sin_3phih_m_phis + sf.F_UTT_sin_phih_m_phis - sf.F_UT_sin_phih_p_phis);
	Real H_121 = -kin.C_1*(sf.F_UT_sin_3phih_m_phis + sf.F_UT_sin_phih_p_phis);
	Real H_222 = kin.C_1*(sf.F_UT_sin_phih_p_phis + sf.F_UTT_sin_phih_m_phis - sf.F_UT_sin_3phih_m_phis);
	H_12 = -H_222;
	H_22 = 4./(sq(kin.lambda_Y)*kin.ph_t_sq)*(
		- kin.lambda_Y*kin.ph_t_sq*kin.Q_sq*H_002
		- sq(kin.lambda_3)*sq(kin.S_x)*H_112
		+ kin.lambda_2*kin.lambda_Y*H_222
		+ 2.*kin.S_x*kin.lambda_3*kin.ph_t*kin.Q*kin.lambda_Y_sqrt*H_012);
	H_32 = 1./kin.ph_t_sq*(H_222 - H_112);
	H_42 = 2./(kin.lambda_Y*kin.ph_t_sq)*(
		kin.lambda_3*kin.S_x*(H_112 - H_222)
		- kin.ph_t*kin.Q*kin.lambda_Y_sqrt*H_012);
	H_61 = 4./(std::pow(kin.lambda_Y_sqrt, 3)*kin.ph_t_sq)*(
		kin.Q*kin.ph_t*kin.lambda_Y_sqrt*H_021
		- kin.lambda_3*kin.S_x*H_121);
	H_81 = 2./(kin.lambda_Y_sqrt*kin.ph_t_sq)*H_121;
}
HadLU::HadLU(Kinematics kin, SfLU sf) {
	Real H_01 = -kin.C_1*sf.F_LU_sin_phih;
	H_50 = (2.*kin.Q)/(kin.ph_t*kin.lambda_Y_sqrt)*H_01;
}
HadLL::HadLL(Kinematics kin, SfLL sf) {
	Real H_023 = -kin.C_1*sf.F_LL_cos_phih;
	Real H_123 = kin.C_1*sf.F_LL;
	H_73 = 4./(std::pow(kin.lambda_Y_sqrt, 3)*kin.ph_t_sq)*(
		kin.Q*kin.ph_t*kin.lambda_Y_sqrt*H_023
		- kin.lambda_3*kin.S_x*H_123);
	H_93 = 2./(kin.lambda_Y_sqrt*kin.ph_t_sq)*H_123;
}
HadLT::HadLT(Kinematics kin, SfLT sf) {
	Real H_012 = -kin.C_1*(sf.F_LT_cos_phis - sf.F_LT_cos_2phih_m_phis);
	Real H_021 = -kin.C_1*(sf.F_LT_cos_2phih_m_phis + sf.F_LT_cos_phis);
	Real H_121 = kin.C_1*sf.F_LT_cos_phih_m_phis;
	H_52 = -(2.*kin.Q)/(kin.ph_t*kin.lambda_Y_sqrt)*H_012;
	H_71 = 4./(std::pow(kin.lambda_Y_sqrt, 3)*kin.ph_t_sq)*(
		kin.Q*kin.ph_t*kin.lambda_Y_sqrt*H_021
		- kin.lambda_3*kin.S_x*H_121);
	H_91 = 2./(kin.lambda_Y_sqrt*kin.ph_t_sq)*H_121;
}
HadUU::HadUU(Kinematics kin, SfSet const& sf) : HadUU(
		kin,
		sf.sf_uu(kin.hadron, kin.x, kin.z, kin.Q_sq, kin.ph_t_sq)) {
	if (sf.target != kin.target) {
		throw TargetMismatch(kin.target, sf.target);
	}
}
HadUL::HadUL(Kinematics kin, SfSet const& sf) : HadUL(
		kin,
		sf.sf_ul(kin.hadron, kin.x, kin.z, kin.Q_sq, kin.ph_t_sq)) {
	if (sf.target != kin.target) {
		throw TargetMismatch(kin.target, sf.target);
	}
}
HadUT::HadUT(Kinematics kin, SfSet const& sf) : HadUT(
		kin,
		sf.sf_ut(kin.hadron, kin.x, kin.z, kin.Q_sq, kin.ph_t_sq)) {
	if (sf.target != kin.target) {
		throw TargetMismatch(kin.target, sf.target);
	}
}
HadLU::HadLU(Kinematics kin, SfSet const& sf) : HadLU(
		kin,
		sf.sf_lu(kin.hadron, kin.x, kin.z, kin.Q_sq, kin.ph_t_sq)) {
	if (sf.target != kin.target) {
		throw TargetMismatch(kin.target, sf.target);
	}
}
HadLL::HadLL(Kinematics kin, SfSet const& sf) : HadLL(
		kin,
		sf.sf_ll(kin.hadron, kin.x, kin.z, kin.Q_sq, kin.ph_t_sq)) {
	if (sf.target != kin.target) {
		throw TargetMismatch(kin.target, sf.target);
	}
}
HadLT::HadLT(Kinematics kin, SfSet const& sf) : HadLT(
		kin,
		sf.sf_lt(kin.hadron, kin.x, kin.z, kin.Q_sq, kin.ph_t_sq)) {
	if (sf.target != kin.target) {
		throw TargetMismatch(kin.target, sf.target);
	}
}

// Radiative hadronic coefficients. Constructed similarly to the standard
// hadronic coefficients, except that the shifted kinematic variables are used
// instead. Also, since the shifted version of the polarization vector is
// involved, all polarized components get mixed into each other, meaning that we
// can't compute, for instance, `HadRadUL`, but must use `HadRadUP`.
HadRadUU::HadRadUU(KinematicsRad kin, SfUU shift_sf) {
	Real H_00 = kin.shift_C_1*shift_sf.F_UUL;
	Real H_01 = -kin.shift_C_1*shift_sf.F_UU_cos_phih;
	Real H_11 = kin.shift_C_1*(shift_sf.F_UU_cos_2phih + shift_sf.F_UUT);
	Real H_22 = kin.shift_C_1*(shift_sf.F_UUT - shift_sf.F_UU_cos_2phih);
	H_10 = H_22;
	H_20 = 4./(sq(kin.shift_lambda_Y)*kin.shift_ph_t_sq)*(
		kin.shift_lambda_Y*kin.shift_ph_t_sq*kin.shift_Q_sq*H_00
		+ sq(kin.shift_lambda_3)*sq(kin.shift_S_x)*H_11
		- kin.shift_lambda_2*kin.shift_lambda_Y*H_22
		-2.*kin.shift_S_x*kin.shift_lambda_3*kin.shift_ph_t*kin.shift_Q*kin.shift_lambda_Y_sqrt*H_01);
	H_30 = 1./kin.shift_ph_t_sq*(H_11 - H_22);
	H_40 = 2./(kin.shift_lambda_Y*kin.shift_ph_t_sq)*(
		kin.shift_lambda_3*kin.shift_S_x*(H_22 - H_11)
		+ kin.shift_ph_t*kin.shift_Q*kin.shift_lambda_Y_sqrt*H_01);
}
HadRadUP::HadRadUP(KinematicsRad kin, SfUP shift_sf) {
	// TODO: Can we fix `shift_rot` being needlessly calculated twice in both
	// `UP` and `LP` versions of `HadRad`?
	Transform3 shift_rot = frame::hadron_from_shift(kin);
	Real H_002 = kin.shift_C_1*shift_sf.F_UTL_sin_phih_m_phis;
	Real H_012 = kin.shift_C_1*(shift_sf.F_UT_sin_phis - shift_sf.F_UT_sin_2phih_m_phis);
	Real H_021 = kin.shift_C_1*(shift_sf.F_UT_sin_2phih_m_phis + shift_sf.F_UT_sin_phis);
	Real H_023 = kin.shift_C_1*shift_sf.F_UL_sin_phih;
	Real H_112 = kin.shift_C_1*(shift_sf.F_UT_sin_3phih_m_phis + shift_sf.F_UTT_sin_phih_m_phis - shift_sf.F_UT_sin_phih_p_phis);
	Real H_121 = -kin.shift_C_1*(shift_sf.F_UT_sin_3phih_m_phis + shift_sf.F_UT_sin_phih_p_phis);
	Real H_123 = -kin.shift_C_1*shift_sf.F_UL_sin_2phih;
	Real H_222 = kin.shift_C_1*(shift_sf.F_UT_sin_phih_p_phis + shift_sf.F_UTT_sin_phih_m_phis - shift_sf.F_UT_sin_3phih_m_phis);
	H_1 = dot(shift_rot, Vec3(
		0.,
		-H_222,
		0.));
	H_2 = dot(shift_rot, Vec3(
		0.,
		4./(sq(kin.shift_lambda_Y)*kin.shift_ph_t_sq)*(
			- kin.shift_lambda_Y*kin.shift_ph_t_sq*kin.shift_Q_sq*H_002
			- sq(kin.shift_lambda_3)*sq(kin.shift_S_x)*H_112
			+ kin.shift_lambda_2*kin.shift_lambda_Y*H_222
			+ 2.*kin.shift_S_x*kin.shift_lambda_3*kin.shift_ph_t*kin.shift_Q*kin.shift_lambda_Y_sqrt*H_012),
		0.));
	H_3 = dot(shift_rot, Vec3(
		0.,
		1./kin.shift_ph_t_sq*(H_222 - H_112),
		0.));
	H_4 = dot(shift_rot, Vec3(
		0.,
		2./(kin.shift_lambda_Y*kin.shift_ph_t_sq)*(
			kin.shift_lambda_3*kin.shift_S_x*(H_112 - H_222)
			- kin.shift_ph_t*kin.shift_Q*kin.shift_lambda_Y_sqrt*H_012),
		0.));
	H_6 = dot(shift_rot, Vec3(
		4./(std::pow(kin.shift_lambda_Y_sqrt, 3)*kin.shift_ph_t_sq)*(
			kin.shift_Q*kin.shift_ph_t*kin.shift_lambda_Y_sqrt*H_021
			- kin.shift_lambda_3*kin.shift_S_x*H_121),
		0.,
		4./(std::pow(kin.shift_lambda_Y_sqrt, 3)*kin.shift_ph_t_sq)*(
			kin.shift_Q*kin.shift_ph_t*kin.shift_lambda_Y_sqrt*H_023
			- kin.shift_lambda_3*kin.shift_S_x*H_123)));
	H_8 = dot(shift_rot, Vec3(
		2./(kin.shift_lambda_Y_sqrt*kin.shift_ph_t_sq)*H_121,
		0.,
		2./(kin.shift_lambda_Y_sqrt*kin.shift_ph_t_sq)*H_123));
}
HadRadLU::HadRadLU(KinematicsRad kin, SfLU shift_sf) {
	Real H_01 = -kin.shift_C_1*shift_sf.F_LU_sin_phih;
	H_50 = (2.*kin.shift_Q)/(kin.shift_ph_t*kin.shift_lambda_Y_sqrt)*H_01;
}
HadRadLP::HadRadLP(KinematicsRad kin, SfLP shift_sf) {
	Transform3 shift_rot = frame::hadron_from_shift(kin);
	Real H_012 = -kin.shift_C_1*(shift_sf.F_LT_cos_phis - shift_sf.F_LT_cos_2phih_m_phis);
	Real H_021 = -kin.shift_C_1*(shift_sf.F_LT_cos_2phih_m_phis + shift_sf.F_LT_cos_phis);
	Real H_023 = -kin.shift_C_1*shift_sf.F_LL_cos_phih;
	Real H_121 = kin.shift_C_1*shift_sf.F_LT_cos_phih_m_phis;
	Real H_123 = kin.shift_C_1*shift_sf.F_LL;
	H_5 = dot(shift_rot, Vec3(
		0.,
		-(2.*kin.shift_Q)/(kin.shift_ph_t*kin.shift_lambda_Y_sqrt)*H_012,
		0.));
	H_7 = dot(shift_rot, Vec3(
		4./(std::pow(kin.shift_lambda_Y_sqrt, 3)*kin.shift_ph_t_sq)*(
			kin.shift_Q*kin.shift_ph_t*kin.shift_lambda_Y_sqrt*H_021
			- kin.shift_lambda_3*kin.shift_S_x*H_121),
		0.,
		4./(std::pow(kin.shift_lambda_Y_sqrt, 3)*kin.shift_ph_t_sq)*(
			kin.shift_Q*kin.shift_ph_t*kin.shift_lambda_Y_sqrt*H_023
			- kin.shift_lambda_3*kin.shift_S_x*H_123)));
	H_9 = dot(shift_rot, Vec3(
		2./(kin.shift_lambda_Y_sqrt*kin.shift_ph_t_sq)*H_121,
		0.,
		2./(kin.shift_lambda_Y_sqrt*kin.shift_ph_t_sq)*H_123));
}
HadRadUU::HadRadUU(KinematicsRad kin, SfSet const& sf) : HadRadUU(
		kin,
		sf.sf_uu(kin.hadron, kin.shift_x, kin.shift_z, kin.shift_Q_sq, kin.shift_ph_t_sq)) {
	if (sf.target != kin.target) {
		throw TargetMismatch(kin.target, sf.target);
	}
}
HadRadUP::HadRadUP(KinematicsRad kin, SfSet const& sf) : HadRadUP(
		kin,
		sf.sf_up(kin.hadron, kin.shift_x, kin.shift_z, kin.shift_Q_sq, kin.shift_ph_t_sq)) {
	if (sf.target != kin.target) {
		throw TargetMismatch(kin.target, sf.target);
	}
}
HadRadLU::HadRadLU(KinematicsRad kin, SfSet const& sf) : HadRadLU(
		kin,
		sf.sf_lu(kin.hadron, kin.shift_x, kin.shift_z, kin.shift_Q_sq, kin.shift_ph_t_sq)) {
	if (sf.target != kin.target) {
		throw TargetMismatch(kin.target, sf.target);
	}
}
HadRadLP::HadRadLP(KinematicsRad kin, SfSet const& sf) : HadRadLP(
		kin,
		sf.sf_lp(kin.hadron, kin.shift_x, kin.shift_z, kin.shift_Q_sq, kin.shift_ph_t_sq)) {
	if (sf.target != kin.target) {
		throw TargetMismatch(kin.target, sf.target);
	}
}

// Infrared-divergent-free radiative hadronic coefficients. Similar to the
// regular hadronic coefficients, except that the difference between the shifted
// and unshifted structure functions is stored. Naively, this difference has an
// issue with catastrophic cancellation for small `R`, but by using the
// appropriate constructors, the infrared-divergent-free radiative coefficients
// can be accurately computed using the diffative of the full radiative
// coefficients near `R=0`.
HadRadFUU::HadRadFUU(KinematicsRad kin, SfUU sf, SfUU shift_sf) {
	HadRadUU had(kin, shift_sf);
	HadUU had_0(kin.project(), sf);
	H_10 = had.H_10;
	H_20 = had.H_20;
	H_30 = had.H_30;
	H_40 = had.H_40;
	H_10_diff = (had.H_10 - had_0.H_10)/kin.R;
	H_20_diff = (had.H_20 - had_0.H_20)/kin.R;
	H_30_diff = (had.H_30 - had_0.H_30)/kin.R;
	H_40_diff = (had.H_40 - had_0.H_40)/kin.R;
}
HadRadFUP::HadRadFUP(KinematicsRad kin, SfUP sf, SfUP shift_sf) {
	HadRadUP had(kin, shift_sf);
	HadUP had_0(kin.project(), sf);
	H_1 = had.H_1;
	H_2 = had.H_2;
	H_3 = had.H_3;
	H_4 = had.H_4;
	H_6 = had.H_6;
	H_8 = had.H_8;
	H_1_diff = (had.H_1 - Vec3(0., had_0.H_12, 0.))/kin.R;
	H_2_diff = (had.H_2 - Vec3(0., had_0.H_22, 0.))/kin.R;
	H_3_diff = (had.H_3 - Vec3(0., had_0.H_32, 0.))/kin.R;
	H_4_diff = (had.H_4 - Vec3(0., had_0.H_42, 0.))/kin.R;
	H_6_diff = (had.H_6 - Vec3(had_0.H_61, 0., had_0.H_63))/kin.R;
	H_8_diff = (had.H_8 - Vec3(had_0.H_81, 0., had_0.H_83))/kin.R;
}
HadRadFLU::HadRadFLU(KinematicsRad kin, SfLU sf, SfLU shift_sf) {
	HadRadLU had(kin, shift_sf);
	HadLU had_0(kin.project(), sf);
	H_50 = had.H_50;
	H_50_diff = (had.H_50 - had_0.H_50)/kin.R;
}
HadRadFLP::HadRadFLP(KinematicsRad kin, SfLP sf, SfLP shift_sf) {
	HadRadLP had(kin, shift_sf);
	HadLP had_0(kin.project(), sf);
	H_5 = had.H_5;
	H_7 = had.H_7;
	H_9 = had.H_9;
	H_5_diff = (had.H_5 - Vec3(0., had_0.H_52, 0.))/kin.R;
	H_7_diff = (had.H_7 - Vec3(had_0.H_71, 0., had_0.H_73))/kin.R;
	H_9_diff = (had.H_9 - Vec3(had_0.H_91, 0., had_0.H_93))/kin.R;
}
HadRadFUU::HadRadFUU(KinematicsRad kin, SfSet const& sf, HadUU had_0) {
	if (sf.target != kin.target) {
		throw TargetMismatch(kin.target, sf.target);
	}
	HadRadUU had(kin, sf);
	H_10 = had.H_10;
	H_20 = had.H_20;
	H_30 = had.H_30;
	H_40 = had.H_40;
	Real delta_R = kin.R_max*DELTA_R_REL;
	if (kin.R < delta_R) {
		Real R_1 = delta_R;
		Real R_2 = 2.*delta_R;
		KinematicsRad kin_1(kin.project(), kin.tau, kin.phi_k, R_1);
		KinematicsRad kin_2(kin.project(), kin.tau, kin.phi_k, R_2);
		HadRadUU had_1(kin_1, sf);
		HadRadUU had_2(kin_2, sf);
		H_10_diff = diff_R(kin.R, delta_R, had_0.H_10, had_1.H_10, had_2.H_10);
		H_20_diff = diff_R(kin.R, delta_R, had_0.H_20, had_1.H_20, had_2.H_20);
		H_30_diff = diff_R(kin.R, delta_R, had_0.H_30, had_1.H_30, had_2.H_30);
		H_40_diff = diff_R(kin.R, delta_R, had_0.H_40, had_1.H_40, had_2.H_40);
	} else {
		H_10_diff = (had.H_10 - had_0.H_10)/kin.R;
		H_20_diff = (had.H_20 - had_0.H_20)/kin.R;
		H_30_diff = (had.H_30 - had_0.H_30)/kin.R;
		H_40_diff = (had.H_40 - had_0.H_40)/kin.R;
	}
}
HadRadFUP::HadRadFUP(KinematicsRad kin, SfSet const& sf, HadUP had_0) {
	if (sf.target != kin.target) {
		throw TargetMismatch(kin.target, sf.target);
	}
	HadRadUP had(kin, sf);
	H_1 = had.H_1;
	H_2 = had.H_2;
	H_3 = had.H_3;
	H_4 = had.H_4;
	H_6 = had.H_6;
	H_8 = had.H_8;
	Vec3 had_0_H_1 = Vec3(0., had_0.H_12, 0.);
	Vec3 had_0_H_2 = Vec3(0., had_0.H_22, 0.);
	Vec3 had_0_H_3 = Vec3(0., had_0.H_32, 0.);
	Vec3 had_0_H_4 = Vec3(0., had_0.H_42, 0.);
	Vec3 had_0_H_6 = Vec3(had_0.H_61, 0., had_0.H_63);
	Vec3 had_0_H_8 = Vec3(had_0.H_81, 0., had_0.H_83);
	Real delta_R = kin.R_max*DELTA_R_REL;
	if (kin.R < delta_R) {
		Real R_1 = delta_R;
		Real R_2 = 2.*delta_R;
		KinematicsRad kin_1(kin.project(), kin.tau, kin.phi_k, R_1);
		KinematicsRad kin_2(kin.project(), kin.tau, kin.phi_k, R_2);
		HadRadUP had_1(kin_1, sf);
		HadRadUP had_2(kin_2, sf);
		H_1_diff = diff_R(kin.R, delta_R, had_0_H_1, had_1.H_1, had_2.H_1);
		H_2_diff = diff_R(kin.R, delta_R, had_0_H_2, had_1.H_2, had_2.H_2);
		H_3_diff = diff_R(kin.R, delta_R, had_0_H_3, had_1.H_3, had_2.H_3);
		H_4_diff = diff_R(kin.R, delta_R, had_0_H_4, had_1.H_4, had_2.H_4);
		H_6_diff = diff_R(kin.R, delta_R, had_0_H_6, had_1.H_6, had_2.H_6);
		H_8_diff = diff_R(kin.R, delta_R, had_0_H_8, had_1.H_8, had_2.H_8);
	} else {
		H_1_diff = (had.H_1 - had_0_H_1)/kin.R;
		H_2_diff = (had.H_2 - had_0_H_2)/kin.R;
		H_3_diff = (had.H_3 - had_0_H_3)/kin.R;
		H_4_diff = (had.H_4 - had_0_H_4)/kin.R;
		H_6_diff = (had.H_6 - had_0_H_6)/kin.R;
		H_8_diff = (had.H_8 - had_0_H_8)/kin.R;
	}
}
HadRadFLU::HadRadFLU(KinematicsRad kin, SfSet const& sf, HadLU had_0) {
	if (sf.target != kin.target) {
		throw TargetMismatch(kin.target, sf.target);
	}
	HadRadLU had(kin, sf);
	H_50 = had.H_50;
	Real delta_R = kin.R_max*DELTA_R_REL;
	if (kin.R < delta_R) {
		Real R_1 = delta_R;
		Real R_2 = 2.*delta_R;
		KinematicsRad kin_1(kin.project(), kin.tau, kin.phi_k, R_1);
		KinematicsRad kin_2(kin.project(), kin.tau, kin.phi_k, R_2);
		HadRadLU had_1(kin_1, sf);
		HadRadLU had_2(kin_2, sf);
		H_50_diff = diff_R(kin.R, delta_R, had_0.H_50, had_1.H_50, had_2.H_50);
	} else {
		H_50_diff = (had.H_50 - had_0.H_50)/kin.R;
	}
}
HadRadFLP::HadRadFLP(KinematicsRad kin, SfSet const& sf, HadLP had_0) {
	if (sf.target != kin.target) {
		throw TargetMismatch(kin.target, sf.target);
	}
	HadRadLP had(kin, sf);
	H_5 = had.H_5;
	H_7 = had.H_7;
	H_9 = had.H_9;
	Vec3 had_0_H_5 = Vec3(0., had_0.H_52, 0.);
	Vec3 had_0_H_7 = Vec3(had_0.H_71, 0., had_0.H_73);
	Vec3 had_0_H_9 = Vec3(had_0.H_91, 0., had_0.H_93);
	Real delta_R = kin.R_max*DELTA_R_REL;
	if (kin.R < delta_R) {
		Real R_1 = delta_R;
		Real R_2 = 2.*delta_R;
		KinematicsRad kin_1(kin.project(), kin.tau, kin.phi_k, R_1);
		KinematicsRad kin_2(kin.project(), kin.tau, kin.phi_k, R_2);
		HadRadLP had_1(kin_1, sf);
		HadRadLP had_2(kin_2, sf);
		H_5_diff = diff_R(kin.R, delta_R, had_0_H_5, had_1.H_5, had_2.H_5);
		H_7_diff = diff_R(kin.R, delta_R, had_0_H_7, had_1.H_7, had_2.H_7);
		H_9_diff = diff_R(kin.R, delta_R, had_0_H_9, had_1.H_9, had_2.H_9);
	} else {
		H_5_diff = (had.H_5 - had_0_H_5)/kin.R;
		H_7_diff = (had.H_7 - had_0_H_7)/kin.R;
		H_9_diff = (had.H_9 - had_0_H_9)/kin.R;
	}
}
HadRadFUU::HadRadFUU(KinematicsRad kin, SfSet const& sf) :
	HadRadFUU(kin, sf, HadUU(kin.project(), sf)) { }
HadRadFUP::HadRadFUP(KinematicsRad kin, SfSet const& sf) :
	HadRadFUP(kin, sf, HadUP(kin.project(), sf)) { }
HadRadFLU::HadRadFLU(KinematicsRad kin, SfSet const& sf) :
	HadRadFLU(kin, sf, HadLU(kin.project(), sf)) { }
HadRadFLP::HadRadFLP(KinematicsRad kin, SfSet const& sf) :
	HadRadFLP(kin, sf, HadLP(kin.project(), sf)) { }

