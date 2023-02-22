#include "sidis/hadronic_coeff.hpp"

#include <cmath>
#include <limits>
#include <stdexcept>
#include <string>

#include "sidis/cut.hpp"
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

// Standard hadronic coefficients. Equation [1.14].
HadBaseUU make_had_base_uu(Kinematics const& kin, SfBaseUU const& sf) {
	HadBaseUU had;
	Real H_00 = kin.C_1*sf.F_UUL;
	Real H_01 = -kin.C_1*sf.F_UU_cos_phih;
	Real H_11 = kin.C_1*(sf.F_UU_cos_2phih + sf.F_UUT);
	Real H_22 = kin.C_1*(sf.F_UUT - sf.F_UU_cos_2phih);
	had.H_10 = H_22;
	had.H_20 = 4./(sq(kin.lambda_Y)*kin.ph_t_sq)*(
		kin.lambda_Y*kin.ph_t_sq*kin.Q_sq*H_00
		+ sq(kin.lambda_3)*sq(kin.S_x)*H_11
		- kin.lambda_2*kin.lambda_Y*H_22
		-2.*kin.S_x*kin.lambda_3*kin.ph_t*kin.Q*kin.lambda_Y_sqrt*H_01);
	had.H_30 = 1./kin.ph_t_sq*(H_11 - H_22);
	had.H_40 = 2./(kin.lambda_Y*kin.ph_t_sq)*(
		kin.lambda_3*kin.S_x*(H_22 - H_11)
		+ kin.ph_t*kin.Q*kin.lambda_Y_sqrt*H_01);
	return had;
}
HadBaseUL make_had_base_ul(Kinematics const& kin, SfBaseUL const& sf) {
	HadBaseUL had;
	Real H_023 = -kin.C_1*sf.F_UL_sin_phih;
	Real H_123 = kin.C_1*sf.F_UL_sin_2phih;
	had.H_63 = 4./(std::pow(kin.lambda_Y_sqrt, 3)*kin.ph_t_sq)*(
		kin.lambda_3*kin.S_x*H_123
		-kin.Q*kin.ph_t*kin.lambda_Y_sqrt*H_023);
	had.H_83 = -2./(kin.lambda_Y_sqrt*kin.ph_t_sq)*H_123;
	return had;
}
HadBaseUT make_had_base_ut(Kinematics const& kin, SfBaseUT const& sf) {
	HadBaseUT had;
	Real H_002 = kin.C_1*sf.F_UTL_sin_phih_m_phis;
	Real H_012 = kin.C_1*(sf.F_UT_sin_phis - sf.F_UT_sin_2phih_m_phis);
	Real H_021 = -kin.C_1*(sf.F_UT_sin_2phih_m_phis + sf.F_UT_sin_phis);
	Real H_112 = kin.C_1*(sf.F_UT_sin_3phih_m_phis + sf.F_UTT_sin_phih_m_phis - sf.F_UT_sin_phih_p_phis);
	Real H_121 = kin.C_1*(sf.F_UT_sin_3phih_m_phis + sf.F_UT_sin_phih_p_phis);
	Real H_222 = kin.C_1*(sf.F_UT_sin_phih_p_phis + sf.F_UTT_sin_phih_m_phis - sf.F_UT_sin_3phih_m_phis);
	had.H_12 = -H_222;
	had.H_22 = 4./(sq(kin.lambda_Y)*kin.ph_t_sq)*(
		- kin.lambda_Y*kin.ph_t_sq*kin.Q_sq*H_002
		- sq(kin.lambda_3)*sq(kin.S_x)*H_112
		+ kin.lambda_2*kin.lambda_Y*H_222
		+ 2.*kin.S_x*kin.lambda_3*kin.ph_t*kin.Q*kin.lambda_Y_sqrt*H_012);
	had.H_32 = 1./kin.ph_t_sq*(H_222 - H_112);
	had.H_42 = 2./(kin.lambda_Y*kin.ph_t_sq)*(
		kin.lambda_3*kin.S_x*(H_112 - H_222)
		- kin.ph_t*kin.Q*kin.lambda_Y_sqrt*H_012);
	had.H_61 = 4./(std::pow(kin.lambda_Y_sqrt, 3)*kin.ph_t_sq)*(
		kin.lambda_3*kin.S_x*H_121
		- kin.Q*kin.ph_t*kin.lambda_Y_sqrt*H_021);
	had.H_81 = -2./(kin.lambda_Y_sqrt*kin.ph_t_sq)*H_121;
	return had;
}
HadBaseUP make_had_base_up(Kinematics const& kin, SfBaseUP const& sf) {
	HadBaseUP had;
	had.ul = make_had_base_ul(kin, sf.ul);
	had.ut = make_had_base_ut(kin, sf.ut);
	return had;
}
HadBaseLU make_had_base_lu(Kinematics const& kin, SfBaseLU const& sf) {
	HadBaseLU had;
	Real H_01 = -kin.C_1*sf.F_LU_sin_phih;
	had.H_50 = (2.*kin.Q)/(kin.ph_t*kin.lambda_Y_sqrt)*H_01;
	return had;
}
HadBaseLL make_had_base_ll(Kinematics const& kin, SfBaseLL const& sf) {
	HadBaseLL had;
	Real H_023 = kin.C_1*sf.F_LL_cos_phih;
	Real H_123 = -kin.C_1*sf.F_LL;
	had.H_73 = 4./(std::pow(kin.lambda_Y_sqrt, 3)*kin.ph_t_sq)*(
		kin.lambda_3*kin.S_x*H_123
		- kin.Q*kin.ph_t*kin.lambda_Y_sqrt*H_023);
	had.H_93 = -2./(kin.lambda_Y_sqrt*kin.ph_t_sq)*H_123;
	return had;
}
HadBaseLT make_had_base_lt(Kinematics const& kin, SfBaseLT const& sf) {
	HadBaseLT had;
	Real H_012 = -kin.C_1*(sf.F_LT_cos_phis - sf.F_LT_cos_2phih_m_phis);
	Real H_021 = kin.C_1*(sf.F_LT_cos_2phih_m_phis + sf.F_LT_cos_phis);
	Real H_121 = -kin.C_1*sf.F_LT_cos_phih_m_phis;
	had.H_52 = -(2.*kin.Q)/(kin.ph_t*kin.lambda_Y_sqrt)*H_012;
	had.H_71 = 4./(std::pow(kin.lambda_Y_sqrt, 3)*kin.ph_t_sq)*(
		kin.lambda_3*kin.S_x*H_121
		- kin.Q*kin.ph_t*kin.lambda_Y_sqrt*H_021);
	had.H_91 = -2./(kin.lambda_Y_sqrt*kin.ph_t_sq)*H_121;
	return had;
}
HadBaseLP make_had_base_lp(Kinematics const& kin, SfBaseLP const& sf) {
	HadBaseLP had;
	had.ll = make_had_base_ll(kin, sf.ll);
	had.lt = make_had_base_lt(kin, sf.lt);
	return had;
}

// Radiative hadronic coefficients. Constructed similarly to the standard
// hadronic coefficients, except that the shifted kinematic variables are used
// instead. Also, since the shifted version of the polarization vector is
// involved, all polarized components get mixed into each other, meaning that we
// can't compute, for instance, `HadRadUL`, but must use `HadRadUP`.
HadRadBaseUU make_had_rad_base_uu(KinematicsRad const& kin, SfBaseUU const& sf) {
	HadRadBaseUU had;
	Real H_00 = kin.shift_C_1*sf.F_UUL;
	Real H_01 = -kin.shift_C_1*sf.F_UU_cos_phih;
	Real H_11 = kin.shift_C_1*(sf.F_UU_cos_2phih + sf.F_UUT);
	Real H_22 = kin.shift_C_1*(sf.F_UUT - sf.F_UU_cos_2phih);
	had.H_10 = H_22;
	had.H_20 = 4./(sq(kin.shift_lambda_Y)*kin.shift_ph_t_sq)*(
		kin.shift_lambda_Y*kin.shift_ph_t_sq*kin.shift_Q_sq*H_00
		+ sq(kin.shift_lambda_3)*sq(kin.shift_S_x)*H_11
		- kin.shift_lambda_2*kin.shift_lambda_Y*H_22
		-2.*kin.shift_S_x*kin.shift_lambda_3*kin.shift_ph_t*kin.shift_Q*kin.shift_lambda_Y_sqrt*H_01);
	had.H_30 = 1./kin.shift_ph_t_sq*(H_11 - H_22);
	had.H_40 = 2./(kin.shift_lambda_Y*kin.shift_ph_t_sq)*(
		kin.shift_lambda_3*kin.shift_S_x*(H_22 - H_11)
		+ kin.shift_ph_t*kin.shift_Q*kin.shift_lambda_Y_sqrt*H_01);
	return had;
}
HadRadBaseUP make_had_rad_base_up(KinematicsRad const& kin, SfBaseUP const& sf) {
	HadRadBaseUP had;
	Transform3 shift_rot = frame::hadron_from_shift(kin);
	Real H_002 = kin.shift_C_1*sf.ut.F_UTL_sin_phih_m_phis;
	Real H_012 = kin.shift_C_1*(sf.ut.F_UT_sin_phis - sf.ut.F_UT_sin_2phih_m_phis);
	Real H_021 = -kin.shift_C_1*(sf.ut.F_UT_sin_2phih_m_phis + sf.ut.F_UT_sin_phis);
	Real H_023 = -kin.shift_C_1*sf.ul.F_UL_sin_phih;
	Real H_112 = kin.shift_C_1*(sf.ut.F_UT_sin_3phih_m_phis + sf.ut.F_UTT_sin_phih_m_phis - sf.ut.F_UT_sin_phih_p_phis);
	Real H_121 = kin.shift_C_1*(sf.ut.F_UT_sin_3phih_m_phis + sf.ut.F_UT_sin_phih_p_phis);
	Real H_123 = kin.shift_C_1*sf.ul.F_UL_sin_2phih;
	Real H_222 = kin.shift_C_1*(sf.ut.F_UT_sin_phih_p_phis + sf.ut.F_UTT_sin_phih_m_phis - sf.ut.F_UT_sin_3phih_m_phis);
	had.H_1 = dot(shift_rot, Vec3(
		0.,
		-H_222,
		0.));
	had.H_2 = dot(shift_rot, Vec3(
		0.,
		4./(sq(kin.shift_lambda_Y)*kin.shift_ph_t_sq)*(
			- kin.shift_lambda_Y*kin.shift_ph_t_sq*kin.shift_Q_sq*H_002
			- sq(kin.shift_lambda_3)*sq(kin.shift_S_x)*H_112
			+ kin.shift_lambda_2*kin.shift_lambda_Y*H_222
			+ 2.*kin.shift_S_x*kin.shift_lambda_3*kin.shift_ph_t*kin.shift_Q*kin.shift_lambda_Y_sqrt*H_012),
		0.));
	had.H_3 = dot(shift_rot, Vec3(
		0.,
		1./kin.shift_ph_t_sq*(H_222 - H_112),
		0.));
	had.H_4 = dot(shift_rot, Vec3(
		0.,
		2./(kin.shift_lambda_Y*kin.shift_ph_t_sq)*(
			kin.shift_lambda_3*kin.shift_S_x*(H_112 - H_222)
			- kin.shift_ph_t*kin.shift_Q*kin.shift_lambda_Y_sqrt*H_012),
		0.));
	had.H_6 = dot(shift_rot, Vec3(
		4./(std::pow(kin.shift_lambda_Y_sqrt, 3)*kin.shift_ph_t_sq)*(
			kin.shift_lambda_3*kin.shift_S_x*H_121
			- kin.shift_Q*kin.shift_ph_t*kin.shift_lambda_Y_sqrt*H_021),
		0.,
		4./(std::pow(kin.shift_lambda_Y_sqrt, 3)*kin.shift_ph_t_sq)*(
			kin.shift_lambda_3*kin.shift_S_x*H_123
			- kin.shift_Q*kin.shift_ph_t*kin.shift_lambda_Y_sqrt*H_023)));
	had.H_8 = dot(shift_rot, Vec3(
		-2./(kin.shift_lambda_Y_sqrt*kin.shift_ph_t_sq)*H_121,
		0.,
		-2./(kin.shift_lambda_Y_sqrt*kin.shift_ph_t_sq)*H_123));
	return had;
}
HadRadBaseLU make_had_rad_base_lu(KinematicsRad const& kin, SfBaseLU const& sf) {
	HadRadBaseLU had;
	Real H_01 = -kin.shift_C_1*sf.F_LU_sin_phih;
	had.H_50 = (2.*kin.shift_Q)/(kin.shift_ph_t*kin.shift_lambda_Y_sqrt)*H_01;
	return had;
}
HadRadBaseLP make_had_rad_base_lp(KinematicsRad const& kin, SfBaseLP const& sf) {
	HadRadBaseLP had;
	Transform3 shift_rot = frame::hadron_from_shift(kin);
	Real H_012 = -kin.shift_C_1*(sf.lt.F_LT_cos_phis - sf.lt.F_LT_cos_2phih_m_phis);
	Real H_021 = kin.shift_C_1*(sf.lt.F_LT_cos_2phih_m_phis + sf.lt.F_LT_cos_phis);
	Real H_023 = kin.shift_C_1*sf.ll.F_LL_cos_phih;
	Real H_121 = -kin.shift_C_1*sf.lt.F_LT_cos_phih_m_phis;
	Real H_123 = -kin.shift_C_1*sf.ll.F_LL;
	had.H_5 = dot(shift_rot, Vec3(
		0.,
		-(2.*kin.shift_Q)/(kin.shift_ph_t*kin.shift_lambda_Y_sqrt)*H_012,
		0.));
	had.H_7 = dot(shift_rot, Vec3(
		4./(std::pow(kin.shift_lambda_Y_sqrt, 3)*kin.shift_ph_t_sq)*(
			kin.shift_lambda_3*kin.shift_S_x*H_121
			- kin.shift_Q*kin.shift_ph_t*kin.shift_lambda_Y_sqrt*H_021),
		0.,
		4./(std::pow(kin.shift_lambda_Y_sqrt, 3)*kin.shift_ph_t_sq)*(
			kin.shift_lambda_3*kin.shift_S_x*H_123
			- kin.shift_Q*kin.shift_ph_t*kin.shift_lambda_Y_sqrt*H_023)));
	had.H_9 = dot(shift_rot, Vec3(
		-2./(kin.shift_lambda_Y_sqrt*kin.shift_ph_t_sq)*H_121,
		0.,
		-2./(kin.shift_lambda_Y_sqrt*kin.shift_ph_t_sq)*H_123));
	return had;
}

// Infrared-divergent-free radiative hadronic coefficients. Similar to the
// regular hadronic coefficients, except that the difference between the shifted
// and unshifted structure functions is stored. Naively, this difference has an
// issue with catastrophic cancellation for small `R`, but by using the
// appropriate constructors, the infrared-divergent-free radiative coefficients
// can be accurately computed using the derivative of the full radiative
// coefficients near `R=0`.
Real const DELTA_R_REL = std::sqrt(2. * std::numeric_limits<Real>::epsilon());
Real diff_R(Real R, Real delta_R, Real f_0, Real f_1, Real f_2) {
	Real R_rel = R/delta_R;
	return 1./(2.*delta_R)*((R_rel - 3.)*f_0 - 2.*(R_rel - 2.)*f_1 + (R_rel - 1.)*f_2);
}
Vec3 diff_R(Real R, Real delta_R, Vec3 f_0, Vec3 f_1, Vec3 f_2) {
	Real R_rel = R/delta_R;
	return 1./(2.*delta_R)*((R_rel - 3.)*f_0 - 2.*(R_rel - 2.)*f_1 + (R_rel - 1.)*f_2);
}

HadRadFBaseUU make_had_rad_f_base_uu_diff(
		Real R,
		HadRadBaseUU const& had, HadBaseUU const& had_0) {
	HadRadFBaseUU had_f;
	had_f.H_10 = had.H_10;
	had_f.H_20 = had.H_20;
	had_f.H_30 = had.H_30;
	had_f.H_40 = had.H_40;
	had_f.H_10_diff = (had.H_10 - had_0.H_10)/R;
	had_f.H_20_diff = (had.H_20 - had_0.H_20)/R;
	had_f.H_30_diff = (had.H_30 - had_0.H_30)/R;
	had_f.H_40_diff = (had.H_40 - had_0.H_40)/R;
	return had_f;
}
HadRadFBaseUU make_had_rad_f_base_uu_deriv(
		Real R, Real delta_R,
		HadRadBaseUU const& had,
		HadBaseUU const& had_0, HadRadBaseUU const& had_1, HadRadBaseUU const& had_2) {
	HadRadFBaseUU had_f;
	had_f.H_10 = had.H_10;
	had_f.H_20 = had.H_20;
	had_f.H_30 = had.H_30;
	had_f.H_40 = had.H_40;
	had_f.H_10_diff = diff_R(R, delta_R, had_0.H_10, had_1.H_10, had_2.H_10);
	had_f.H_20_diff = diff_R(R, delta_R, had_0.H_20, had_1.H_20, had_2.H_20);
	had_f.H_30_diff = diff_R(R, delta_R, had_0.H_30, had_1.H_30, had_2.H_30);
	had_f.H_40_diff = diff_R(R, delta_R, had_0.H_40, had_1.H_40, had_2.H_40);
	return had_f;
}
HadRadFBaseUP make_had_rad_f_base_up_diff(
		Real R,
		HadRadBaseUP const& had, HadBaseUP const& had_0) {
	HadRadFBaseUP had_f;
	had_f.H_1 = had.H_1;
	had_f.H_2 = had.H_2;
	had_f.H_3 = had.H_3;
	had_f.H_4 = had.H_4;
	had_f.H_6 = had.H_6;
	had_f.H_8 = had.H_8;
	had_f.H_1_diff = (had.H_1 - Vec3(0., had_0.ut.H_12, 0.))/R;
	had_f.H_2_diff = (had.H_2 - Vec3(0., had_0.ut.H_22, 0.))/R;
	had_f.H_3_diff = (had.H_3 - Vec3(0., had_0.ut.H_32, 0.))/R;
	had_f.H_4_diff = (had.H_4 - Vec3(0., had_0.ut.H_42, 0.))/R;
	had_f.H_6_diff = (had.H_6 - Vec3(had_0.ut.H_61, 0., had_0.ul.H_63))/R;
	had_f.H_8_diff = (had.H_8 - Vec3(had_0.ut.H_81, 0., had_0.ul.H_83))/R;
	return had_f;
}
HadRadFBaseUP make_had_rad_f_base_up_deriv(
		Real R, Real delta_R,
		HadRadBaseUP const& had,
		HadBaseUP const& had_0, HadRadBaseUP const& had_1, HadRadBaseUP const& had_2) {
	HadRadFBaseUP had_f;
	had_f.H_1 = had.H_1;
	had_f.H_2 = had.H_2;
	had_f.H_3 = had.H_3;
	had_f.H_4 = had.H_4;
	had_f.H_6 = had.H_6;
	had_f.H_8 = had.H_8;
	had_f.H_1_diff = diff_R(R, delta_R, Vec3(0., had_0.ut.H_12, 0.), had_1.H_1, had_2.H_1);
	had_f.H_2_diff = diff_R(R, delta_R, Vec3(0., had_0.ut.H_22, 0.), had_1.H_2, had_2.H_2);
	had_f.H_3_diff = diff_R(R, delta_R, Vec3(0., had_0.ut.H_32, 0.), had_1.H_3, had_2.H_3);
	had_f.H_4_diff = diff_R(R, delta_R, Vec3(0., had_0.ut.H_42, 0.), had_1.H_4, had_2.H_4);
	had_f.H_6_diff = diff_R(R, delta_R, Vec3(had_0.ut.H_61, 0., had_0.ul.H_63), had_1.H_6, had_2.H_6);
	had_f.H_8_diff = diff_R(R, delta_R, Vec3(had_0.ut.H_81, 0., had_0.ul.H_83), had_1.H_8, had_2.H_8);
	return had_f;
}
HadRadFBaseLU make_had_rad_f_base_lu_diff(
		Real R,
		HadRadBaseLU const& had, HadBaseLU const& had_0) {
	HadRadFBaseLU had_f;
	had_f.H_50 = had.H_50;
	had_f.H_50_diff = (had.H_50 - had_0.H_50)/R;
	return had_f;
}
HadRadFBaseLU make_had_rad_f_base_lu_deriv(
		Real R, Real delta_R,
		HadRadBaseLU const& had,
		HadBaseLU const& had_0, HadRadBaseLU const& had_1, HadRadBaseLU const& had_2) {
	HadRadFBaseLU had_f;
	had_f.H_50 = had.H_50;
	had_f.H_50_diff = diff_R(R, delta_R, had_0.H_50, had_1.H_50, had_2.H_50);
	return had_f;
}
HadRadFBaseLP make_had_rad_f_base_lp_diff(
		Real R,
		HadRadBaseLP const& had, HadBaseLP const& had_0) {
	HadRadFBaseLP had_f;
	had_f.H_5 = had.H_5;
	had_f.H_7 = had.H_7;
	had_f.H_9 = had.H_9;
	had_f.H_5_diff = (had.H_5 - Vec3(0., had_0.lt.H_52, 0.))/R;
	had_f.H_7_diff = (had.H_7 - Vec3(had_0.lt.H_71, 0., had_0.ll.H_73))/R;
	had_f.H_9_diff = (had.H_9 - Vec3(had_0.lt.H_91, 0., had_0.ll.H_93))/R;
	return had_f;
}
HadRadFBaseLP make_had_rad_f_base_lp_deriv(
		Real R, Real delta_R,
		HadRadBaseLP const& had,
		HadBaseLP const& had_0, HadRadBaseLP const& had_1, HadRadBaseLP const& had_2) {
	HadRadFBaseLP had_f;
	had_f.H_5 = had.H_5;
	had_f.H_7 = had.H_7;
	had_f.H_9 = had.H_9;
	had_f.H_5_diff = diff_R(R, delta_R, Vec3(0., had_0.lt.H_52, 0.), had_1.H_5, had_2.H_5);
	had_f.H_7_diff = diff_R(R, delta_R, Vec3(had_0.lt.H_71, 0., had_0.ll.H_73), had_1.H_7, had_2.H_7);
	had_f.H_9_diff = diff_R(R, delta_R, Vec3(had_0.lt.H_91, 0., had_0.ll.H_93), had_1.H_9, had_2.H_9);
	return had_f;
}

}

// Below here, the various constructors for the hadronic coefficients are
// defined using the above functions. The constructors try to avoid making more
// calls than necessary to SfSet that provides the structure functions.
// TODO: Verify at run-time that the kinematics has the same target as the
// structure function set.

// HadBaseXX constructors.
HadBaseUU::HadBaseUU(Kinematics const& kin, SfSet const& sf_set) :
	HadBaseUU(make_had_base_uu(
		kin,
		sf_set.sf_base_uu(kin.hadron, kin.x, kin.z, kin.Q_sq, kin.ph_t_sq))) { }
HadBaseUL::HadBaseUL(Kinematics const& kin, SfSet const& sf_set) :
	HadBaseUL(make_had_base_ul(
		kin,
		sf_set.sf_base_ul(kin.hadron, kin.x, kin.z, kin.Q_sq, kin.ph_t_sq))) { }
HadBaseUT::HadBaseUT(Kinematics const& kin, SfSet const& sf_set) :
	HadBaseUT(make_had_base_ut(
		kin,
		sf_set.sf_base_ut(kin.hadron, kin.x, kin.z, kin.Q_sq, kin.ph_t_sq))) { }
HadBaseUP::HadBaseUP(Kinematics const& kin, SfSet const& sf_set) :
	HadBaseUP(make_had_base_up(
		kin,
		sf_set.sf_base_up(kin.hadron, kin.x, kin.z, kin.Q_sq, kin.ph_t_sq))) { }
HadBaseLU::HadBaseLU(Kinematics const& kin, SfSet const& sf_set) :
	HadBaseLU(make_had_base_lu(
		kin,
		sf_set.sf_base_lu(kin.hadron, kin.x, kin.z, kin.Q_sq, kin.ph_t_sq))) { }
HadBaseLL::HadBaseLL(Kinematics const& kin, SfSet const& sf_set) :
	HadBaseLL(make_had_base_ll(
		kin,
		sf_set.sf_base_ll(kin.hadron, kin.x, kin.z, kin.Q_sq, kin.ph_t_sq))) { }
HadBaseLT::HadBaseLT(Kinematics const& kin, SfSet const& sf_set) :
	HadBaseLT(make_had_base_lt(
		kin,
		sf_set.sf_base_lt(kin.hadron, kin.x, kin.z, kin.Q_sq, kin.ph_t_sq))) { }
HadBaseLP::HadBaseLP(Kinematics const& kin, SfSet const& sf_set) :
	HadBaseLP(make_had_base_lp(
		kin,
		sf_set.sf_base_lp(kin.hadron, kin.x, kin.z, kin.Q_sq, kin.ph_t_sq))) { }

// HadXX constructors.
HadUU::HadUU(Kinematics const& kin, SfSet const& sf_set) {
	SfUU sf = sf_set.sf_uu(kin.hadron, kin.x, kin.z, kin.Q_sq, kin.ph_t_sq);
	uu = make_had_base_uu(kin, sf.uu);
}
HadUL::HadUL(Kinematics const& kin, SfSet const& sf_set) {
	SfUL sf = sf_set.sf_ul(kin.hadron, kin.x, kin.z, kin.Q_sq, kin.ph_t_sq);
	uu = make_had_base_uu(kin, sf.uu);
	ul = make_had_base_ul(kin, sf.ul);
}
HadUT::HadUT(Kinematics const& kin, SfSet const& sf_set) {
	SfUT sf = sf_set.sf_ut(kin.hadron, kin.x, kin.z, kin.Q_sq, kin.ph_t_sq);
	uu = make_had_base_uu(kin, sf.uu);
	ut = make_had_base_ut(kin, sf.ut);
}
HadUP::HadUP(Kinematics const& kin, SfSet const& sf_set) {
	SfUP sf = sf_set.sf_up(kin.hadron, kin.x, kin.z, kin.Q_sq, kin.ph_t_sq);
	uu = make_had_base_uu(kin, sf.uu);
	ul = make_had_base_ul(kin, sf.ul);
	ut = make_had_base_ut(kin, sf.ut);
}
HadLU::HadLU(Kinematics const& kin, SfSet const& sf_set) {
	SfLU sf = sf_set.sf_lu(kin.hadron, kin.x, kin.z, kin.Q_sq, kin.ph_t_sq);
	uu = make_had_base_uu(kin, sf.uu);
	lu = make_had_base_lu(kin, sf.lu);
}
HadLL::HadLL(Kinematics const& kin, SfSet const& sf_set) {
	SfLL sf = sf_set.sf_ll(kin.hadron, kin.x, kin.z, kin.Q_sq, kin.ph_t_sq);
	uu = make_had_base_uu(kin, sf.uu);
	ul = make_had_base_ul(kin, sf.ul);
	lu = make_had_base_lu(kin, sf.lu);
	ll = make_had_base_ll(kin, sf.ll);
}
HadLT::HadLT(Kinematics const& kin, SfSet const& sf_set) {
	SfLT sf = sf_set.sf_lt(kin.hadron, kin.x, kin.z, kin.Q_sq, kin.ph_t_sq);
	uu = make_had_base_uu(kin, sf.uu);
	ut = make_had_base_ut(kin, sf.ut);
	lu = make_had_base_lu(kin, sf.lu);
	lt = make_had_base_lt(kin, sf.lt);
}
HadLP::HadLP(Kinematics const& kin, SfSet const& sf_set) {
	SfLP sf = sf_set.sf_lp(kin.hadron, kin.x, kin.z, kin.Q_sq, kin.ph_t_sq);
	uu = make_had_base_uu(kin, sf.uu);
	ul = make_had_base_ul(kin, sf.ul);
	ut = make_had_base_ut(kin, sf.ut);
	lu = make_had_base_lu(kin, sf.lu);
	ll = make_had_base_ll(kin, sf.ll);
	lt = make_had_base_lt(kin, sf.lt);
}

// HadRadBaseXX constructors.
HadRadBaseUU::HadRadBaseUU(KinematicsRad const& kin, SfSet const& sf_set) :
	HadRadBaseUU(make_had_rad_base_uu(
		kin,
		sf_set.sf_base_uu(kin.hadron, kin.shift_x, kin.shift_z, kin.shift_Q_sq, kin.shift_ph_t_sq))) { }
HadRadBaseUP::HadRadBaseUP(KinematicsRad const& kin, SfSet const& sf_set) :
	HadRadBaseUP(make_had_rad_base_up(
		kin,
		sf_set.sf_base_up(kin.hadron, kin.shift_x, kin.shift_z, kin.shift_Q_sq, kin.shift_ph_t_sq))) { }
HadRadBaseLU::HadRadBaseLU(KinematicsRad const& kin, SfSet const& sf_set) :
	HadRadBaseLU(make_had_rad_base_lu(
		kin,
		sf_set.sf_base_lu(kin.hadron, kin.shift_x, kin.shift_z, kin.shift_Q_sq, kin.shift_ph_t_sq))) { }
HadRadBaseLP::HadRadBaseLP(KinematicsRad const& kin, SfSet const& sf_set) :
	HadRadBaseLP(make_had_rad_base_lp(
		kin,
		sf_set.sf_base_lp(kin.hadron, kin.shift_x, kin.shift_z, kin.shift_Q_sq, kin.shift_ph_t_sq))) { }

// HadRadXX constructors.
HadRadUU::HadRadUU(KinematicsRad const& kin, SfSet const& sf_set) {
	SfUU sf = sf_set.sf_uu(kin.hadron, kin.shift_x, kin.shift_z, kin.shift_Q_sq, kin.shift_ph_t_sq);
	uu = make_had_rad_base_uu(kin, sf.uu);
}
HadRadUP::HadRadUP(KinematicsRad const& kin, SfSet const& sf_set) {
	SfUP sf = sf_set.sf_up(kin.hadron, kin.shift_x, kin.shift_z, kin.shift_Q_sq, kin.shift_ph_t_sq);
	uu = make_had_rad_base_uu(kin, sf.uu);
	up = make_had_rad_base_up(kin, { sf.ul, sf.ut });
}
HadRadLU::HadRadLU(KinematicsRad const& kin, SfSet const& sf_set) {
	SfLU sf = sf_set.sf_lu(kin.hadron, kin.shift_x, kin.shift_z, kin.shift_Q_sq, kin.shift_ph_t_sq);
	uu = make_had_rad_base_uu(kin, sf.uu);
	lu = make_had_rad_base_lu(kin, sf.lu);
}
HadRadLP::HadRadLP(KinematicsRad const& kin, SfSet const& sf_set) {
	SfLP sf = sf_set.sf_lp(kin.hadron, kin.shift_x, kin.shift_z, kin.shift_Q_sq, kin.shift_ph_t_sq);
	uu = make_had_rad_base_uu(kin, sf.uu);
	up = make_had_rad_base_up(kin, { sf.ul, sf.ut });
	lu = make_had_rad_base_lu(kin, sf.lu);
	lp = make_had_rad_base_lp(kin, { sf.ll, sf.lt });
}

// HadRadFBaseXX constructors.
HadRadFBaseUU::HadRadFBaseUU(KinematicsRad const& kin, SfSet const& sf_set, HadBaseUU const& had_0) {
	Kinematics kin_0 = kin.project();
	HadRadBaseUU had(kin, sf_set);
	Real delta_R = cut::R_bound(kin_0, kin.tau, kin.phi_k).max()*DELTA_R_REL;
	if (kin.R < delta_R) {
		KinematicsRad kin_1(kin_0, kin.tau, kin.phi_k, delta_R);
		KinematicsRad kin_2(kin_0, kin.tau, kin.phi_k, 2.*delta_R);
		HadRadBaseUU had_1(kin_1, sf_set);
		HadRadBaseUU had_2(kin_2, sf_set);
		*this = make_had_rad_f_base_uu_deriv(kin.R, delta_R, had, had_0, had_1, had_2);
	} else {
		*this = make_had_rad_f_base_uu_diff(kin.R, had, had_0);
	}
}
HadRadFBaseUP::HadRadFBaseUP(KinematicsRad const& kin, SfSet const& sf_set, HadBaseUP const& had_0) {
	Kinematics kin_0 = kin.project();
	HadRadBaseUP had(kin, sf_set);
	Real delta_R = cut::R_bound(kin_0, kin.tau, kin.phi_k).max()*DELTA_R_REL;
	if (kin.R < delta_R) {
		KinematicsRad kin_1(kin_0, kin.tau, kin.phi_k, delta_R);
		KinematicsRad kin_2(kin_0, kin.tau, kin.phi_k, 2.*delta_R);
		HadRadBaseUP had_1(kin_1, sf_set);
		HadRadBaseUP had_2(kin_2, sf_set);
		*this = make_had_rad_f_base_up_deriv(kin.R, delta_R, had, had_0, had_1, had_2);
	} else {
		*this = make_had_rad_f_base_up_diff(kin.R, had, had_0);
	}
}
HadRadFBaseLU::HadRadFBaseLU(KinematicsRad const& kin, SfSet const& sf_set, HadBaseLU const& had_0) {
	Kinematics kin_0 = kin.project();
	HadRadBaseLU had(kin, sf_set);
	Real delta_R = cut::R_bound(kin_0, kin.tau, kin.phi_k).max()*DELTA_R_REL;
	if (kin.R < delta_R) {
		KinematicsRad kin_1(kin_0, kin.tau, kin.phi_k, delta_R);
		KinematicsRad kin_2(kin_0, kin.tau, kin.phi_k, 2.*delta_R);
		HadRadBaseLU had_1(kin_1, sf_set);
		HadRadBaseLU had_2(kin_2, sf_set);
		*this = make_had_rad_f_base_lu_deriv(kin.R, delta_R, had, had_0, had_1, had_2);
	} else {
		*this = make_had_rad_f_base_lu_diff(kin.R, had, had_0);
	}
}
HadRadFBaseLP::HadRadFBaseLP(KinematicsRad const& kin, SfSet const& sf_set, HadBaseLP const& had_0) {
	Kinematics kin_0 = kin.project();
	HadRadBaseLP had(kin, sf_set);
	Real delta_R = cut::R_bound(kin_0, kin.tau, kin.phi_k).max()*DELTA_R_REL;
	if (kin.R < delta_R) {
		KinematicsRad kin_1(kin_0, kin.tau, kin.phi_k, delta_R);
		KinematicsRad kin_2(kin_0, kin.tau, kin.phi_k, 2.*delta_R);
		HadRadBaseLP had_1(kin_1, sf_set);
		HadRadBaseLP had_2(kin_2, sf_set);
		*this = make_had_rad_f_base_lp_deriv(kin.R, delta_R, had, had_0, had_1, had_2);
	} else {
		*this = make_had_rad_f_base_lp_diff(kin.R, had, had_0);
	}
}

// HadRadFXX constructors.
HadRadFUU::HadRadFUU(KinematicsRad const& kin, SfSet const& sf_set, HadUU const& had_0) {
	Kinematics kin_0 = kin.project();
	HadRadUU had(kin, sf_set);
	Real delta_R = cut::R_bound(kin_0, kin.tau, kin.phi_k).max()*DELTA_R_REL;
	if (kin.R < delta_R) {
		KinematicsRad kin_1(kin_0, kin.tau, kin.phi_k, delta_R);
		KinematicsRad kin_2(kin_0, kin.tau, kin.phi_k, 2.*delta_R);
		HadRadUU had_1(kin_1, sf_set);
		HadRadUU had_2(kin_2, sf_set);
		uu = make_had_rad_f_base_uu_deriv(kin.R, delta_R, had.uu, had_0.uu, had_1.uu, had_2.uu);
	} else {
		uu = make_had_rad_f_base_uu_diff(kin.R, had.uu, had_0.uu);
	}
}
HadRadFUP::HadRadFUP(KinematicsRad const& kin, SfSet const& sf_set, HadUP const& had_0) {
	Kinematics kin_0 = kin.project();
	HadBaseUP had_0_up(had_0.ul, had_0.ut);
	HadRadUP had(kin, sf_set);
	Real delta_R = cut::R_bound(kin_0, kin.tau, kin.phi_k).max()*DELTA_R_REL;
	if (kin.R < delta_R) {
		KinematicsRad kin_1(kin_0, kin.tau, kin.phi_k, delta_R);
		KinematicsRad kin_2(kin_0, kin.tau, kin.phi_k, 2.*delta_R);
		HadRadUP had_1(kin_1, sf_set);
		HadRadUP had_2(kin_2, sf_set);
		uu = make_had_rad_f_base_uu_deriv(kin.R, delta_R, had.uu, had_0.uu, had_1.uu, had_2.uu);
		up = make_had_rad_f_base_up_deriv(kin.R, delta_R, had.up, had_0_up, had_1.up, had_2.up);
	} else {
		uu = make_had_rad_f_base_uu_diff(kin.R, had.uu, had_0.uu);
		up = make_had_rad_f_base_up_diff(kin.R, had.up, had_0_up);
	}
}
HadRadFLU::HadRadFLU(KinematicsRad const& kin, SfSet const& sf_set, HadLU const& had_0) {
	Kinematics kin_0 = kin.project();
	HadRadLU had(kin, sf_set);
	Real delta_R = cut::R_bound(kin_0, kin.tau, kin.phi_k).max()*DELTA_R_REL;
	if (kin.R < delta_R) {
		KinematicsRad kin_1(kin_0, kin.tau, kin.phi_k, delta_R);
		KinematicsRad kin_2(kin_0, kin.tau, kin.phi_k, 2.*delta_R);
		HadRadLU had_1(kin_1, sf_set);
		HadRadLU had_2(kin_2, sf_set);
		uu = make_had_rad_f_base_uu_deriv(kin.R, delta_R, had.uu, had_0.uu, had_1.uu, had_2.uu);
		lu = make_had_rad_f_base_lu_deriv(kin.R, delta_R, had.lu, had_0.lu, had_1.lu, had_2.lu);
	} else {
		uu = make_had_rad_f_base_uu_diff(kin.R, had.uu, had_0.uu);
		lu = make_had_rad_f_base_lu_diff(kin.R, had.lu, had_0.lu);
	}
}
HadRadFLP::HadRadFLP(KinematicsRad const& kin, SfSet const& sf_set, HadLP const& had_0) {
	Kinematics kin_0 = kin.project();
	HadBaseUP had_0_up(had_0.ul, had_0.ut);
	HadBaseLP had_0_lp(had_0.ll, had_0.lt);
	HadRadLP had(kin, sf_set);
	Real delta_R = cut::R_bound(kin_0, kin.tau, kin.phi_k).max()*DELTA_R_REL;
	if (kin.R < delta_R) {
		KinematicsRad kin_1(kin_0, kin.tau, kin.phi_k, delta_R);
		KinematicsRad kin_2(kin_0, kin.tau, kin.phi_k, 2.*delta_R);
		HadRadLP had_1(kin_1, sf_set);
		HadRadLP had_2(kin_2, sf_set);
		uu = make_had_rad_f_base_uu_deriv(kin.R, delta_R, had.uu, had_0.uu, had_1.uu, had_2.uu);
		up = make_had_rad_f_base_up_deriv(kin.R, delta_R, had.up, had_0_up, had_1.up, had_2.up);
		lu = make_had_rad_f_base_lu_deriv(kin.R, delta_R, had.lu, had_0.lu, had_1.lu, had_2.lu);
		lp = make_had_rad_f_base_lp_deriv(kin.R, delta_R, had.lp, had_0_lp, had_1.lp, had_2.lp);
	} else {
		uu = make_had_rad_f_base_uu_diff(kin.R, had.uu, had_0.uu);
		up = make_had_rad_f_base_up_diff(kin.R, had.up, had_0_up);
		lu = make_had_rad_f_base_lu_diff(kin.R, had.lu, had_0.lu);
		lp = make_had_rad_f_base_lp_diff(kin.R, had.lp, had_0_lp);
	}
}

