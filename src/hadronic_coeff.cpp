#include "sidis/hadronic_coeff.hpp"

#include <cmath>

#include "sidis/kinematics.hpp"
#include "sidis/structure_function.hpp"
#include "sidis/extra/math.hpp"

using namespace sidis;
using namespace sidis::had;
using namespace sidis::kin;
using namespace sidis::math;
using namespace sidis::sf;

HadUU::HadUU(Kinematics kin, SfUU sf) {
	Real H_00 = sf.F_UUL;
	Real H_01 = -sf.F_UU_cos_phih;
	Real H_11 = sf.F_UU_cos_2phih + sf.F_UUT;
	Real H_22 = sf.F_UUT - sf.F_UU_cos_2phih;
	H_1 = H_22;
	H_2 = 4./(sq(kin.lambda_Y)*kin.ph_t_sq)*(
		kin.lambda_Y*kin.ph_t_sq*kin.Q_sq*H_00
		+ sq(kin.lambda_3)*sq(kin.S_x)*H_11
		- kin.lambda_2*kin.lambda_Y*H_22
		-2.*kin.S_x*kin.lambda_3*kin.ph_t*kin.Q*kin.lambda_Y_sqrt*H_01);
	H_3 = 1./kin.ph_t_sq*(H_11 - H_22);
	H_4 = 2./(kin.lambda_Y*kin.ph_t_sq)*(
		kin.lambda_3*kin.S_x*(H_22 - H_11)
		+ kin.ph_t*kin.Q*kin.lambda_Y_sqrt*H_01);
}

HadUL::HadUL(Kinematics kin, SfUL sf) {
	Real H_023 = sf.F_UL_sin_phih;
	Real H_123 = -sf.F_UL_sin_2phih;
	H_6 = 4./(std::pow(kin.lambda_Y_sqrt, 3)*kin.ph_t_sq)*(
		kin.Q*kin.ph_t*kin.lambda_Y_sqrt*H_023
		- kin.lambda_3*kin.S_x*H_123);
	H_8 = 2./(kin.lambda_Y_sqrt*kin.ph_t_sq)*H_123;
}

HadUT1::HadUT1(Kinematics kin, SfUT sf) {
	Real H_021 = sf.F_UT_sin_2phih_m_phis + sf.F_UT_sin_phis;
	Real H_121 = -(sf.F_UT_sin_3phih_m_phis + sf.F_UT_sin_phih_p_phis);
	H_6 = 4./(std::pow(kin.lambda_Y_sqrt, 3)*kin.ph_t_sq)*(
		kin.Q*kin.ph_t*kin.lambda_Y_sqrt*H_021
		- kin.lambda_3*kin.S_x*H_121);
	H_8 = 2./(kin.lambda_Y_sqrt*kin.ph_t_sq)*H_121;
}

HadUT2::HadUT2(Kinematics kin, SfUT sf) {
	Real H_002 = sf.F_UTL_sin_phih_m_phis;
	Real H_012 = sf.F_UT_sin_phis - sf.F_UT_sin_2phih_m_phis;
	Real H_112 = sf.F_UT_sin_3phih_m_phis + sf.F_UTT_sin_phih_m_phis - sf.F_UT_sin_phih_p_phis;
	Real H_222 = sf.F_UT_sin_phih_p_phis + sf.F_UTT_sin_phih_m_phis - sf.F_UT_sin_3phih_m_phis;
	H_1 = -H_222;
	H_2 = 4./(sq(kin.lambda_Y)*kin.ph_t_sq)*(
		- kin.lambda_Y*kin.ph_t_sq*kin.Q_sq*H_002
		- sq(kin.lambda_3)*sq(kin.S_x)*H_112
		+ kin.lambda_2*kin.lambda_Y*H_222
		+ 2.*kin.S_x*kin.lambda_3*kin.ph_t*kin.Q*kin.lambda_Y_sqrt*H_012);
	H_3 = 1./kin.ph_t_sq*(H_222 - H_112);
	H_4 = 2./(kin.lambda_Y*kin.ph_t_sq)*(
		kin.lambda_3*kin.S_x*(H_112 - H_222)
		- kin.ph_t*kin.Q*kin.lambda_Y_sqrt*H_012);
}

HadLU::HadLU(Kinematics kin, SfLU sf) {
	Real H_01 = -sf.F_LU_sin_phih;
	H_5 = (2.*kin.Q)/(kin.ph_t*kin.lambda_Y_sqrt)*H_01;
}

HadLL::HadLL(Kinematics kin, SfLL sf) {
	Real H_023 = -sf.F_LL_cos_phih;
	Real H_123 = sf.F_LL;
	H_7 = 4./(std::pow(kin.lambda_Y_sqrt, 3)*kin.ph_t_sq)*(
		kin.Q*kin.ph_t*kin.lambda_Y_sqrt*H_023
		- kin.lambda_3*kin.S_x*H_123);
	H_9 = 2./(kin.lambda_Y_sqrt*kin.ph_t_sq)*H_123;
}

HadLT1::HadLT1(Kinematics kin, SfLT sf) {
	Real H_021 = -(sf.F_LT_cos_2phih_m_phis + sf.F_LT_cos_phis);
	Real H_121 = sf.F_LT_cos_phih_m_phis;
	H_7 = 4./(std::pow(kin.lambda_Y_sqrt, 3)*kin.ph_t_sq)*(
		kin.Q*kin.ph_t*kin.lambda_Y_sqrt*H_021
		- kin.lambda_3*kin.S_x*H_121);
	H_9 = 2./(kin.lambda_Y_sqrt*kin.ph_t_sq)*H_121;
}

HadLT2::HadLT2(Kinematics kin, SfLT sf) {
	Real H_012 = -(sf.F_LT_cos_phis - sf.F_LT_cos_2phih_m_phis);
	H_5 = -(2.*kin.Q)/(kin.ph_t*kin.lambda_Y_sqrt)*H_012;
}

