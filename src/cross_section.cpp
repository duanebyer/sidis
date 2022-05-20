#include "sidis/cross_section.hpp"

#include <cmath>
#include <limits>
#include <string>

#include "sidis/bound.hpp"
#include "sidis/constant.hpp"
#include "sidis/cut.hpp"
#include "sidis/frame.hpp"
#include "sidis/hadronic_coeff.hpp"
#include "sidis/kinematics.hpp"
#include "sidis/leptonic_coeff.hpp"
#include "sidis/structure_function.hpp"
#include "sidis/vector.hpp"
#include "sidis/extra/integrate.hpp"
#include "sidis/extra/math.hpp"

using namespace sidis;
using namespace sidis::cut;
using namespace sidis::had;
using namespace sidis::kin;
using namespace sidis::lep;
using namespace sidis::math;
using namespace sidis::sf;
using namespace sidis::xs;

// Macro that computes the cross-section from the base cross-sections in an
// optimized way. For example, if the polarization is zero, then the
// cross-section can be computed just using the UU base cross-section.
#define SIDIS_MACRO_XS_FROM_BASE(name, Lep, Had, kin, sf, b, lambda_e, eta) \
	[&](){ \
		Real uu = 0.; \
		Vec3 up = Vec3::ZERO; \
		Real lu = 0.; \
		Vec3 lp = Vec3::ZERO; \
		Lep##UU lep_uu((kin)); \
		Had##UU had_uu((kin), (sf)); \
		uu = name##_uu_base((b), lep_uu, had_uu); \
		if (eta.z != 0. && (eta.x != 0. || eta.y != 0.)) { \
			Lep##UP lep_up((kin)); \
			Had##UL had_ul((kin), (sf)); \
			Had##UT had_ut((kin), (sf)); \
			up.z = name##_ul_base((b), lep_up, had_ul); \
			up.x = eta.x == 0. ? 0. : name##_ut1_base((b), lep_up, had_ut); \
			up.y = eta.y == 0. ? 0. : name##_ut2_base((b), lep_uu, had_ut); \
		} else if (eta.z != 0.) { \
			Lep##UP lep_up((kin)); \
			Had##UL had_ul((kin), (sf)); \
			up.z = name##_ul_base((b), lep_up, had_ul); \
		} else { \
			Had##UT had_ut((kin), (sf)); \
			if (eta.x != 0.) { \
				Lep##UP lep_up((kin)); \
				up.x = name##_ut1_base((b), lep_up, had_ut); \
			} \
			up.y = eta.y == 0. ? 0. : name##_ut2_base((b), lep_uu, had_ut); \
		} \
		if (lambda_e != 0.) { \
			Lep##LU lep_lu((kin)); \
			Had##LU had_lu((kin), (sf)); \
			lu = name##_lu_base((b), lep_lu, had_lu); \
			if (eta.z != 0. && (eta.x != 0. || eta.y != 0.)) { \
				Lep##LP lep_lp((kin)); \
				Had##LL had_ll((kin), (sf)); \
				Had##LT had_lt((kin), (sf)); \
				lp.z = name##_ll_base((b), lep_lp, had_ll); \
				lp.x = eta.x == 0. ? 0. : name##_lt1_base((b), lep_lp, had_lt); \
				lp.y = eta.y == 0. ? 0. : name##_lt2_base((b), lep_lu, had_lt); \
			} else if (eta.z != 0.) { \
				Lep##LP lep_lp((kin)); \
				Had##LL had_ll((kin), (sf)); \
				lp.z = name##_ll_base((b), lep_lp, had_ll); \
			} else { \
				Had##LT had_lt((kin), (sf)); \
				if (eta.x != 0.) { \
					Lep##LP lep_lp((kin)); \
					lp.x = name##_lt1_base((b), lep_lp, had_lt); \
				} \
				lp.y = eta.y == 0. ? 0. : name##_lt2_base((b), lep_lu, had_lt); \
			} \
		} \
		return uu + dot(up, eta) + (lambda_e)*(lu + dot(lp, (eta))); \
	}()

// Similar to `SIDIS_MACRO_XS_FROM_BASE`, except this one works with base cross-
// sections where the XL, XT1, and XT2 cases are all grouped together into an
// XP case.
#define SIDIS_MACRO_XS_FROM_BASE_P(name, Lep, Had, kin, sf, b, lambda_e, eta) \
	[&](){ \
		Real uu = 0.; \
		Vec3 up = Vec3::ZERO; \
		Real lu = 0.; \
		Vec3 lp = Vec3::ZERO; \
		Lep##UU lep_uu((kin)); \
		Had##UU had_uu((kin), (sf)); \
		uu = name##_uu_base((b), lep_uu, had_uu); \
		if (eta != Vec3::ZERO) { \
			Lep##UX lep_ux(lep_uu, Lep##UP(kin)); \
			Had##UP had_up((kin), (sf)); \
			up = name##_up_base((b), lep_ux, had_up); \
		} \
		if (lambda_e != 0.) { \
			Lep##LU lep_lu((kin)); \
			Had##LU had_lu((kin), (sf)); \
			lu = name##_lu_base((b), lep_lu, had_lu); \
			if (eta != Vec3::ZERO) { \
				Lep##LX lep_lx(lep_lu, Lep##LP(kin)); \
				Had##LP had_lp((kin), (sf)); \
				lp = name##_lp_base((b), lep_lx, had_lp); \
			} \
		} \
		return uu + dot(up, eta) + (lambda_e)*(lu + dot(lp, (eta))); \
	}()

namespace {

Real delta_vert_rad_0(Kinematics const& kin) {
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

	// Equation [1.40].
	Real rho = 1./lambda_m_sqrt*(
		(Q_m_sq + lambda_m_sqrt)*S_prime
		- 2.*sq(kin.m)*X_prime);
	Real S_phi = Q_m_sq/lambda_m_sqrt*(
		lambda_S_prime*sq(L_S_prime)/4. - lambda_X_prime*sq(L_X_prime)/4.
		+ dilog(1. - 1./(sum_S_prime*S_prime)*rho)
		+ dilog(1. - (sum_S_prime*S_prime)/(4.*sq(kin.m)*kin.mx_sq)*rho)
		- dilog(1. - (sum_X_prime*X_prime)/(kin.mx_sq*sq(sum_m)*kin.Q_sq)*rho)
		- dilog(1. - (4.*sq(kin.m))/(sum_X_prime*X_prime*sq(sum_m)*kin.Q_sq)*rho));

	// Equation [1.52].
	Real delta = 0.5*S_prime*L_S_prime + 0.5*X_prime*L_X_prime + S_phi - 2.
		+ (1.5*kin.Q_sq + 4.*sq(kin.m))*L_m
		- Q_m_sq/lambda_m_sqrt*(
			0.5*lambda_m*sq(L_m)
			+ 2.*dilog((2.*lambda_m_sqrt)/(kin.Q_sq + lambda_m_sqrt))
			- 0.5*sq(PI));
	return delta;
}

}

Real xs::born(Kinematics const& kin, SfSet const& sf, Real lambda_e, Vec3 eta) {
	Born b(kin);
	return SIDIS_MACRO_XS_FROM_BASE(born, LepBorn, Had, kin, sf, b, lambda_e, eta);
}

Real xs::amm(Kinematics const& kin, SfSet const& sf, Real lambda_e, Vec3 eta) {
	Amm b(kin);
	return SIDIS_MACRO_XS_FROM_BASE(amm, LepAmm, Had, kin, sf, b, lambda_e, eta);
}

Real xs::nrad_ir(Kinematics const& kin, SfSet const& sf, Real lambda_e, Vec3 eta, Real k_0_bar) {
	Nrad b(kin, k_0_bar);
	return SIDIS_MACRO_XS_FROM_BASE(nrad_ir, LepNrad, Had, kin, sf, b, lambda_e, eta);
}

Real xs::rad(KinematicsRad const& kin, SfSet const& sf, Real lambda_e, Vec3 eta) {
	Rad b(kin);
	return SIDIS_MACRO_XS_FROM_BASE_P(rad, LepRad, HadRad, kin, sf, b, lambda_e, eta);
}

Real xs::rad_f(KinematicsRad const& kin, SfSet const& sf, Real lambda_e, Vec3 eta) {
	Rad b(kin);
	return SIDIS_MACRO_XS_FROM_BASE_P(rad_f, LepRad, HadRadF, kin, sf, b, lambda_e, eta);
}

EstErr xs::nrad_integ(Kinematics const& kin, SfSet const& sf, Real lambda_e, Vec3 eta, Real k_0_bar, IntegParams params) {
	// The soft part of the radiative cross-section (below `k_0_bar`) is bundled
	// into the return value here.
	Real xs_nrad_ir = nrad_ir(kin, sf, lambda_e, eta, k_0_bar);
	// TODO: The integration parameters should be modified here to account for
	// the `xs_nrad_ir` contribution.
	EstErr xs_rad_f = rad_f_integ(kin, sf, lambda_e, eta, k_0_bar, params);
	return { xs_nrad_ir + xs_rad_f.val, xs_rad_f.err };
}

EstErr xs::rad_f_integ(Kinematics const& kin, SfSet const& sf, Real lambda_e, Vec3 eta, Real k_0_bar, IntegParams params) {
	HadXX had_0(kin, sf);
	CutRad cut;
	cut.k_0_bar = Bound(0., k_0_bar);
	EstErr xs_integ = integrate<3>(
		[&](std::array<Real, 3> x) {
			KinematicsRad kin_rad;
			Real jacobian;
			if (!take(cut, kin, x.data(), &kin_rad, &jacobian)) {
				return 0.;
			}

			Rad b(kin_rad);
			LepRadXX lep(kin_rad);
			HadRadFXX had(kin_rad, sf, had_0);
			Real uu = rad_f_uu_base(b, lep, had);
			Vec3 up = rad_f_up_base(b, lep, had);
			Real lu = rad_f_lu_base(b, lep, had);
			Vec3 lp = rad_f_lp_base(b, lep, had);
			Real xs = uu + dot(eta, up) + lambda_e*(lu + dot(eta, lp));
			if (std::isnan(xs)) {
				// If the result is `NaN`, it most likely means we went out of
				// the allowed region for the structure function grids (or we
				// are in a kinematically disallowed region). In that case, just
				// return zero.
				return 0.;
			} else {
				return jacobian * xs;
			}
		},
		std::array<Real, 3>{ 0., 0., 0. },
		std::array<Real, 3>{ 1., 1., 1. },
		params);
	return xs_integ;
}

EstErr xs::rad_integ(Kinematics const& kin, SfSet const& sf, Real lambda_e, Vec3 eta, Real k_0_bar, IntegParams params) {
	CutRad cut;
	cut.k_0_bar = Bound(k_0_bar, INF);
	EstErr xs_integ = integrate<3>(
		[&](std::array<Real, 3> x) {
			KinematicsRad kin_rad;
			Real jacobian;
			if (!take(cut, kin, x.data(), &kin_rad, &jacobian)) {
				return 0.;
			}

			Rad b(kin_rad);
			LepRadXX lep(kin_rad);
			HadRadXX had(kin_rad, sf);
			Real uu = rad_uu_base(b, lep, had);
			Vec3 up = rad_up_base(b, lep, had);
			Real lu = rad_lu_base(b, lep, had);
			Vec3 lp = rad_lp_base(b, lep, had);
			Real xs = uu + dot(eta, up) + lambda_e*(lu + dot(eta, lp));
			if (std::isnan(xs)) {
				return 0.;
			} else {
				return jacobian * xs;
			}
		},
		std::array<Real, 3>{ 0., 0., 0. },
		std::array<Real, 3>{ 1., 1., 1. },
		params);
	return xs_integ;
}

// Radiative corrections to Born cross-section.
Real xs::delta_vert_rad_ir(Kinematics const& kin, Real k_0_bar) {
	// Paragraph following equation [1.C17].
	Real k0_max = (kin.mx_sq - sq(kin.Mth))/(2.*kin.mx);
	if (!(k_0_bar > 0.)) {
		return -INFINITY;
	}
	Real Q_m_sq = kin.Q_sq + 2.*sq(kin.m);
	Real diff_m = sqrt1p_1m((4.*sq(kin.m))/kin.Q_sq);
	Real sum_m = 2. + diff_m;
	Real lambda_m = kin.Q_sq*(kin.Q_sq + 4.*sq(kin.m));
	Real lambda_m_sqrt = std::sqrt(lambda_m);
	Real L_m = 1./lambda_m_sqrt*std::log(sum_m/diff_m);
	Real delta_0 = delta_vert_rad_0(kin);
	// This comes from subtracting `delta_H` (equation [1.38]) from `delta_VR`
	// (equation [1.52]).
	Real delta_shift = 2.*(Q_m_sq*L_m - 1.)*std::log(
		k_0_bar < k0_max ?
		(2.*k_0_bar)/kin.m :
		(kin.mx_sq - sq(kin.Mth))/(kin.m*kin.mx));
	return delta_0 + delta_shift;
}
Real xs::delta_rad_ir_hard(Kinematics const& kin, Real k_0_bar) {
	Real k0_max = (kin.mx_sq - sq(kin.Mth))/(2.*kin.mx);
	if (!(k_0_bar > 0.)) {
		return INFINITY;
	} else if (!(k_0_bar < k0_max)) {
		return 0.;
	}
	Real Q_m_sq = kin.Q_sq + 2.*sq(kin.m);
	Real diff_m = sqrt1p_1m((4.*sq(kin.m))/kin.Q_sq);
	Real sum_m = 2. + diff_m;
	Real lambda_m = kin.Q_sq*(kin.Q_sq + 4.*sq(kin.m));
	Real lambda_m_sqrt = std::sqrt(lambda_m);
	Real L_m = 1./lambda_m_sqrt*std::log(sum_m/diff_m);
	// Equation [1.38].
	Real delta = 2.*(Q_m_sq*L_m - 1.)*std::log(
		(kin.mx_sq - sq(kin.Mth))/(2.*k_0_bar*kin.mx));
	return delta;
}

Real xs::delta_vac_lep(Kinematics const& kin) {
	// Equation [1.50].
	Real ms[3] = { MASS_E, MASS_MU, MASS_TAU };
	Real delta = 0.;
	for (unsigned idx = 0; idx < 3; ++idx) {
		Real m = ms[idx];
		Real lambda_sqrt = std::sqrt(kin.Q_sq*(kin.Q_sq + 4.*sq(m)));
		Real diff_m = sqrt1p_1m((4.*sq(m))/kin.Q_sq);
		Real sum_m = 2. + diff_m;
		Real L_m = 1./lambda_sqrt*std::log(sum_m/diff_m);
		delta += 2./3.L*(kin.Q_sq
			+ 2.*sq(m))*L_m
			- 10./9.L
			+ (8.*sq(m))/(3.*kin.Q_sq)*(1. - 2.*sq(m)*L_m);
	}
	return delta;
}

Real xs::delta_vac_had(Kinematics const& kin) {
	if (kin.Q_sq < 1.) {
		return -(2.*PI)/ALPHA*(-1.345e-9L - 2.302e-3L*std::log(1. + 4.091L*kin.Q_sq));
	} else if (kin.Q_sq < 64.) {
		return -(2.*PI)/ALPHA*(-1.512e-3L - 2.822e-3L*std::log(1. + 1.218L*kin.Q_sq));
	} else {
		return -(2.*PI)/ALPHA*(-1.1344e-3L - 3.0680e-3L*std::log(1. + 0.99992L*kin.Q_sq));
	}
}

// Born base functions.
Born::Born(Kinematics const& kin) :
	// Equation [1.15]. The `Q^4` factor has been absorbed into `C_1`.
	coeff((sq(ALPHA)*kin.S*sq(kin.S_x))/(8.*kin.M*kin.ph_l*kin.lambda_S)) { }

Real xs::born_uu_base(Born const& b, LepBornUU const& lep, HadUU const& had) {
	return b.coeff*(
		lep.theta_1*had.H_10
		+ lep.theta_2*had.H_20
		+ lep.theta_3*had.H_30
		+ lep.theta_4*had.H_40);
}
Real xs::born_ul_base(Born const& b, LepBornUP const& lep, HadUL const& had) {
	return b.coeff*(lep.theta_6*had.H_63 + lep.theta_8*had.H_83);
}
Real xs::born_ut1_base(Born const& b, LepBornUP const& lep, HadUT const& had) {
	return b.coeff*(lep.theta_6*had.H_61 + lep.theta_8*had.H_81);
}
Real xs::born_ut2_base(Born const& b, LepBornUU const& lep, HadUT const& had) {
	return b.coeff*(
		lep.theta_1*had.H_12
		+ lep.theta_2*had.H_22
		+ lep.theta_3*had.H_32
		+ lep.theta_4*had.H_42);
}
Vec3 xs::born_up_base(Born const& b, LepBornUX const& lep, HadUP const& had) {
	return Vec3(
		born_ut1_base(b, lep, had),
		born_ut2_base(b, lep, had),
		born_ul_base(b, lep, had));
}
Real xs::born_lu_base(Born const& b, LepBornLU const& lep, HadLU const& had) {
	return b.coeff*lep.theta_5*had.H_50;
}
Real xs::born_ll_base(Born const& b, LepBornLP const& lep, HadLL const& had) {
	return b.coeff*(lep.theta_7*had.H_73 + lep.theta_9*had.H_93);
}
Real xs::born_lt1_base(Born const& b, LepBornLP const& lep, HadLT const& had) {
	return b.coeff*(lep.theta_7*had.H_71 + lep.theta_9*had.H_91);
}
Real xs::born_lt2_base(Born const& b, LepBornLU const& lep, HadLT const& had) {
	return b.coeff*lep.theta_5*had.H_52;
}
Vec3 xs::born_lp_base(Born const& b, LepBornLX const& lep, HadLP const& had) {
	return Vec3(
		born_lt1_base(b, lep, had),
		born_lt2_base(b, lep, had),
		born_ll_base(b, lep, had));
}

// AMM base functions.
Amm::Amm(Kinematics const& kin) {
	// Equation [1.53]. The `Q^4` factor has been absorbed into `C_1`.
	Real lambda_m = kin.Q_sq*(kin.Q_sq + 4.*sq(kin.m));
	Real lambda_m_sqrt = std::sqrt(lambda_m);
	Real diff_m = sqrt1p_1m((4.*sq(kin.m))/kin.Q_sq);
	Real sum_m = 2. + diff_m;
	Real L_m = 1./lambda_m_sqrt*std::log(sum_m/diff_m);
	coeff = L_m*kin.Q_sq*(std::pow(ALPHA, 3)*sq(kin.m)*kin.S*sq(kin.S_x))
		/(16.*PI*kin.M*kin.ph_l*kin.lambda_S);
}

Real xs::amm_uu_base(Amm const& b, LepAmmUU const& lep, HadUU const& had) {
	return b.coeff*(
		lep.theta_1*had.H_10
		+ lep.theta_2*had.H_20
		+ lep.theta_3*had.H_30
		+ lep.theta_4*had.H_40);
}
Real xs::amm_ul_base(Amm const& b, LepAmmUP const& lep, HadUL const& had) {
	return b.coeff*(lep.theta_6*had.H_63 + lep.theta_8*had.H_83);
}
Real xs::amm_ut1_base(Amm const& b, LepAmmUP const& lep, HadUT const& had) {
	return b.coeff*(lep.theta_6*had.H_61 + lep.theta_8*had.H_81);
}
Real xs::amm_ut2_base(Amm const& b, LepAmmUU const& lep, HadUT const& had) {
	return b.coeff*(
		lep.theta_1*had.H_12
		+ lep.theta_2*had.H_22
		+ lep.theta_3*had.H_32
		+ lep.theta_4*had.H_42);
}
Vec3 xs::amm_up_base(Amm const& b, LepAmmUX const& lep, HadUP const& had) {
	return Vec3(
		amm_ut1_base(b, lep, had),
		amm_ut2_base(b, lep, had),
		amm_ul_base(b, lep, had));
}
Real xs::amm_lu_base(Amm const& b, LepAmmLU const& lep, HadLU const& had) {
	return b.coeff*lep.theta_5*had.H_50;
}
Real xs::amm_ll_base(Amm const& b, LepAmmLP const& lep, HadLL const& had) {
	return b.coeff*(lep.theta_7*had.H_73 + lep.theta_9*had.H_93);
}
Real xs::amm_lt1_base(Amm const& b, LepAmmLP const& lep, HadLT const& had) {
	return b.coeff*(lep.theta_7*had.H_71 + lep.theta_9*had.H_91);
}
Real xs::amm_lt2_base(Amm const& b, LepAmmLU const& lep, HadLT const& had) {
	return b.coeff*lep.theta_5*had.H_52;
}
Vec3 xs::amm_lp_base(Amm const& b, LepAmmLX const& lep, HadLP const& had) {
	return Vec3(
		amm_lt1_base(b, lep, had),
		amm_lt2_base(b, lep, had),
		amm_ll_base(b, lep, had));
}

// Non-radiative infrared-divergence-free base functions.
Nrad::Nrad(Kinematics const& kin, Real k_0_bar) {
	Born born(kin);
	Amm amm(kin);
	Real born_factor = 1. + ALPHA/PI*(
		delta_vert_rad_ir(kin, k_0_bar)
		+ delta_vac_lep(kin)
		+ delta_vac_had(kin));
	coeff_born = born_factor*born.coeff;
	coeff_amm = amm.coeff;
}

Real xs::nrad_ir_uu_base(Nrad const& b, LepNradUU const& lep, HadUU const& had) {
	return
		(b.coeff_born*lep.born_uu.theta_1 + b.coeff_amm*lep.amm_uu.theta_1)*had.H_10
		+ (b.coeff_born*lep.born_uu.theta_2 + b.coeff_amm*lep.amm_uu.theta_2)*had.H_20
		+ (b.coeff_born*lep.born_uu.theta_3 + b.coeff_amm*lep.amm_uu.theta_3)*had.H_30
		+ (b.coeff_born*lep.born_uu.theta_4 + b.coeff_amm*lep.amm_uu.theta_4)*had.H_40;
}
Real xs::nrad_ir_ul_base(Nrad const& b, LepNradUP const& lep, HadUL const& had) {
	return
		(b.coeff_born*lep.born_up.theta_6 + b.coeff_amm*lep.amm_up.theta_6)*had.H_63
		+ (b.coeff_born*lep.born_up.theta_8 + b.coeff_amm*lep.amm_up.theta_8)*had.H_83;
}
Real xs::nrad_ir_ut1_base(Nrad const& b, LepNradUP const& lep, HadUT const& had) {
	return
		(b.coeff_born*lep.born_up.theta_6 + b.coeff_amm*lep.amm_up.theta_6)*had.H_61
		+ (b.coeff_born*lep.born_up.theta_8 + b.coeff_amm*lep.amm_up.theta_8)*had.H_81;
}
Real xs::nrad_ir_ut2_base(Nrad const& b, LepNradUU const& lep, HadUT const& had) {
	return
		(b.coeff_born*lep.born_uu.theta_1 + b.coeff_amm*lep.amm_uu.theta_1)*had.H_12
		+ (b.coeff_born*lep.born_uu.theta_2 + b.coeff_amm*lep.amm_uu.theta_2)*had.H_22
		+ (b.coeff_born*lep.born_uu.theta_3 + b.coeff_amm*lep.amm_uu.theta_3)*had.H_32
		+ (b.coeff_born*lep.born_uu.theta_4 + b.coeff_amm*lep.amm_uu.theta_4)*had.H_42;
}
Vec3 xs::nrad_ir_up_base(Nrad const& b, LepNradUX const& lep, HadUP const& had) {
	return Vec3(
		nrad_ir_ut1_base(b, lep, had),
		nrad_ir_ut2_base(b, lep, had),
		nrad_ir_ul_base(b, lep, had));
}
Real xs::nrad_ir_lu_base(Nrad const& b, LepNradLU const& lep, HadLU const& had) {
	return (b.coeff_born*lep.born_lu.theta_5 + b.coeff_amm*lep.amm_lu.theta_5)*had.H_50;
}
Real xs::nrad_ir_ll_base(Nrad const& b, LepNradLP const& lep, HadLL const& had) {
	return
		(b.coeff_born*lep.born_lp.theta_7 + b.coeff_amm*lep.amm_lp.theta_7)*had.H_73
		+ (b.coeff_born*lep.born_lp.theta_9 + b.coeff_amm*lep.amm_lp.theta_9)*had.H_93;
}
Real xs::nrad_ir_lt1_base(Nrad const& b, LepNradLP const& lep, HadLT const& had) {
	return
		(b.coeff_born*lep.born_lp.theta_7 + b.coeff_amm*lep.amm_lp.theta_7)*had.H_71
		+ (b.coeff_born*lep.born_lp.theta_9 + b.coeff_amm*lep.amm_lp.theta_9)*had.H_91;
}
Real xs::nrad_ir_lt2_base(Nrad const& b, LepNradLU const& lep, HadLT const& had) {
	return (b.coeff_born*lep.born_lu.theta_5 + b.coeff_amm*lep.amm_lu.theta_5)*had.H_52;
}
Vec3 xs::nrad_ir_lp_base(Nrad const& b, LepNradLX const& lep, HadLP const& had) {
	return Vec3(
		nrad_ir_lt1_base(b, lep, had),
		nrad_ir_lt2_base(b, lep, had),
		nrad_ir_ll_base(b, lep, had));
}

// Radiative base functions.
Rad::Rad(KinematicsRad const& kin) {
	// Equation [1.43].
	coeff = -(std::pow(ALPHA, 3)*kin.S*sq(kin.S_x))
		/(64.*sq(PI)*kin.M*kin.ph_l*kin.lambda_S*kin.lambda_Y_sqrt);
	R = kin.R;
}

Real xs::rad_uu_base(Rad const& b, LepRadUU const& lep, HadRadUU const& had) {
	return b.coeff*(
		1./b.R*(
			lep.theta_011*had.H_10
			+ lep.theta_021*had.H_20
			+ lep.theta_031*had.H_30
			+ lep.theta_041*had.H_40)
		+ (
			lep.theta_012*had.H_10
			+ lep.theta_022*had.H_20
			+ lep.theta_032*had.H_30
			+ lep.theta_042*had.H_40)
		+ b.R*(
			lep.theta_013*had.H_10
			+ lep.theta_023*had.H_20
			+ lep.theta_033*had.H_30
			+ lep.theta_043*had.H_40));
}
Vec3 xs::rad_up_base(Rad const& b, LepRadUX const& lep, HadRadUP const& had) {
	return b.coeff*(
		1./b.R*(
			lep.theta_011*had.H_1
			+ lep.theta_021*had.H_2
			+ lep.theta_031*had.H_3
			+ lep.theta_041*had.H_4
			+ lep.theta_061*had.H_6
			+ lep.theta_081*had.H_8)
		+ (
			lep.theta_012*had.H_1
			+ lep.theta_022*had.H_2
			+ lep.theta_032*had.H_3
			+ lep.theta_042*had.H_4
			+ lep.theta_062*had.H_6
			+ lep.theta_082*had.H_8)
		+ b.R*(
			lep.theta_013*had.H_1
			+ lep.theta_023*had.H_2
			+ lep.theta_033*had.H_3
			+ lep.theta_043*had.H_4
			+ lep.theta_063*had.H_6
			+ lep.theta_083*had.H_8)
		+ b.R*b.R*(
			lep.theta_064*had.H_6
			+ lep.theta_084*had.H_8));
}
Real xs::rad_lu_base(Rad const& b, LepRadLU const& lep, HadRadLU const& had) {
	return b.coeff*(
		1./b.R*(lep.theta_051 + lep.theta_151)*had.H_50
		+ (lep.theta_052 + lep.theta_152)*had.H_50
		+ b.R*(lep.theta_053 + lep.theta_153)*had.H_50);
}
Vec3 xs::rad_lp_base(Rad const& b, LepRadLX const& lep, HadRadLP const& had) {
	return b.coeff*(
		1./b.R*(
			(lep.theta_051 + lep.theta_151)*had.H_5
			+ (lep.theta_071 + lep.theta_171)*had.H_7
			+ (lep.theta_091 + lep.theta_191)*had.H_9)
		+ (
		 	(lep.theta_052 + lep.theta_152)*had.H_5
			+ (lep.theta_072 + lep.theta_172)*had.H_7
			+ (lep.theta_092 + lep.theta_192)*had.H_9)
		+ b.R*(
		 	(lep.theta_053 + lep.theta_153)*had.H_5
			+ (lep.theta_073 + lep.theta_173)*had.H_7
			+ (lep.theta_093 + lep.theta_193)*had.H_9)
		+ b.R*b.R*(
			(lep.theta_074 + lep.theta_174)*had.H_7
			+ (lep.theta_094 + lep.theta_194)*had.H_9));
}

Real xs::rad_f_uu_base(Rad const& b, LepRadUU const& lep, HadRadFUU const& had) {
	return b.coeff*(
		(
			lep.theta_011*had.H_10_diff
			+ lep.theta_021*had.H_20_diff
			+ lep.theta_031*had.H_30_diff
			+ lep.theta_041*had.H_40_diff)
		+ (
			lep.theta_012*had.H_10
			+ lep.theta_022*had.H_20
			+ lep.theta_032*had.H_30
			+ lep.theta_042*had.H_40)
		+ b.R*(
			lep.theta_013*had.H_10
			+ lep.theta_023*had.H_20
			+ lep.theta_033*had.H_30
			+ lep.theta_043*had.H_40));
}
Vec3 xs::rad_f_up_base(Rad const& b, LepRadUX const& lep, HadRadFUP const& had) {
	return b.coeff*(
		(
			lep.theta_011*had.H_1_diff
			+ lep.theta_021*had.H_2_diff
			+ lep.theta_031*had.H_3_diff
			+ lep.theta_041*had.H_4_diff
			+ lep.theta_061*had.H_6_diff
			+ lep.theta_081*had.H_8_diff)
		+ (
			lep.theta_012*had.H_1
			+ lep.theta_022*had.H_2
			+ lep.theta_032*had.H_3
			+ lep.theta_042*had.H_4
			+ lep.theta_062*had.H_6
			+ lep.theta_082*had.H_8)
		+ b.R*(
			lep.theta_013*had.H_1
			+ lep.theta_023*had.H_2
			+ lep.theta_033*had.H_3
			+ lep.theta_043*had.H_4
			+ lep.theta_063*had.H_6
			+ lep.theta_083*had.H_8)
		+ b.R*b.R*(
			lep.theta_064*had.H_6
			+ lep.theta_084*had.H_8));
}
Real xs::rad_f_lu_base(Rad const& b, LepRadLU const& lep, HadRadFLU const& had) {
	return b.coeff*(
		(lep.theta_051 + lep.theta_151)*had.H_50_diff
		+ (lep.theta_052 + lep.theta_152)*had.H_50
		+ b.R*(lep.theta_053 + lep.theta_153)*had.H_50);
}
Vec3 xs::rad_f_lp_base(Rad const& b, LepRadLX const& lep, HadRadFLP const& had) {
	return b.coeff*(
		(
			(lep.theta_051 + lep.theta_151)*had.H_5_diff
			+ (lep.theta_071 + lep.theta_171)*had.H_7_diff
			+ (lep.theta_091 + lep.theta_191)*had.H_9_diff)
		+ (
		 	(lep.theta_052 + lep.theta_152)*had.H_5
			+ (lep.theta_072 + lep.theta_172)*had.H_7
			+ (lep.theta_092 + lep.theta_192)*had.H_9)
		+ b.R*(
		 	(lep.theta_053 + lep.theta_153)*had.H_5
			+ (lep.theta_073 + lep.theta_173)*had.H_7
			+ (lep.theta_093 + lep.theta_193)*had.H_9)
		+ b.R*b.R*(
			(lep.theta_074 + lep.theta_174)*had.H_7
			+ (lep.theta_094 + lep.theta_194)*had.H_9));
}

