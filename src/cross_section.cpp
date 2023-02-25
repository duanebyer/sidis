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
#include "sidis/phenom.hpp"
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
using namespace sidis::ph;
using namespace sidis::sf;
using namespace sidis::xs;

// Macro that computes the cross-section from the base cross-sections in an
// optimized way. For example, if the polarization is zero, then the
// cross-section can be computed just using the UU base cross-section.
#define SIDIS_MACRO_XS_FROM_BASE(name, Lep, Had, kin, sf, b, lambda_e, eta) ([&]() { \
	/* Create a mask describing the polarization state. */ \
	unsigned pol_mask = (((lambda_e) != 0.) << 3) \
		| (((eta).x != 0.) << 2) \
		| (((eta).y != 0.) << 1) \
		| (((eta).z != 0.) << 0); \
	Real uu = 0.; \
	Vec3 up = VEC3_ZERO; \
	Real lu = 0.; \
	Vec3 lp = VEC3_ZERO; \
	switch (pol_mask) { \
	case 0:  /* 0000 */ \
		{ \
			Lep##UU lep((kin)); \
			Had##UU had((kin), (sf)); \
			uu = name##_base_uu((b), lep.uu, had.uu); \
		} \
		break; \
	case 1:  /* 0001 */ \
		{ \
			Lep##UP lep((kin)); \
			Had##UL had((kin), (sf)); \
			uu = name##_base_uu((b), lep.uu, had.uu); \
			up.z = name##_base_ul((b), lep.up, had.ul); \
		} \
		break; \
	case 2:  /* 0010 */ \
		{ \
			Lep##UU lep((kin)); \
			Had##UT had((kin), (sf)); \
			uu = name##_base_uu((b), lep.uu, had.uu); \
			up.y = name##_base_ut2((b), lep.uu, had.ut); \
		} \
		break; \
	case 3:  /* 0011 */ \
		{ \
			Lep##UP lep((kin)); \
			Had##UP had((kin), (sf)); \
			uu = name##_base_uu((b), lep.uu, had.uu); \
			up.y = name##_base_ut2((b), lep.uu, had.ut); \
			up.z = name##_base_ul((b), lep.up, had.ul); \
		} \
		break; \
	case 4:  /* 0100 */ \
		{ \
			Lep##UP lep((kin)); \
			Had##UT had((kin), (sf)); \
			uu = name##_base_uu((b), lep.uu, had.uu); \
			up.x = name##_base_ut1((b), lep.up, had.ut); \
		} \
		break; \
	case 5:  /* 0101 */ \
		{ \
			Lep##UP lep((kin)); \
			Had##UP had((kin), (sf)); \
			uu = name##_base_uu((b), lep.uu, had.uu); \
			up.x = name##_base_ut1((b), lep.up, had.ut); \
			up.z = name##_base_ul((b), lep.up, had.ul); \
		} \
		break; \
	case 6:  /* 0110 */ \
		{ \
			Lep##UP lep((kin)); \
			Had##UT had((kin), (sf)); \
			uu = name##_base_uu((b), lep.uu, had.uu); \
			up.x = name##_base_ut1((b), lep.up, had.ut); \
			up.y = name##_base_ut2((b), lep.uu, had.ut); \
		} \
		break; \
	case 7:  /* 0111 */ \
		{ \
			Lep##UP lep((kin)); \
			Had##UP had((kin), (sf)); \
			uu = name##_base_uu((b), lep.uu, had.uu); \
			up.x = name##_base_ut1((b), lep.up, had.ut); \
			up.y = name##_base_ut2((b), lep.uu, had.ut); \
			up.z = name##_base_ul((b), lep.up, had.ul); \
		} \
		break; \
	case 8:  /* 1000 */ \
		{ \
			Lep##LU lep((kin)); \
			Had##LU had((kin), (sf)); \
			uu = name##_base_uu((b), lep.uu, had.uu); \
			lu = name##_base_lu((b), lep.lu, had.lu); \
		} \
		break; \
	case 9:  /* 1001 */ \
		{ \
			Lep##LP lep((kin)); \
			Had##LL had((kin), (sf)); \
			uu = name##_base_uu((b), lep.uu, had.uu); \
			up.z = name##_base_ul((b), lep.up, had.ul); \
			lu = name##_base_lu((b), lep.lu, had.lu); \
			lp.z = name##_base_ll((b), lep.lp, had.ll); \
		} \
		break; \
	case 10: /* 1010 */ \
		{ \
			Lep##LU lep((kin)); \
			Had##LT had((kin), (sf)); \
			uu = name##_base_uu((b), lep.uu, had.uu); \
			up.y = name##_base_ut2((b), lep.uu, had.ut); \
			lu = name##_base_lu((b), lep.lu, had.lu); \
			lp.y = name##_base_lt2((b), lep.lu, had.lt); \
		} \
		break; \
	case 11: /* 1011 */ \
		{ \
			Lep##LP lep((kin)); \
			Had##LP had((kin), (sf)); \
			uu = name##_base_uu((b), lep.uu, had.uu); \
			up.y = name##_base_ut2((b), lep.uu, had.ut); \
			up.z = name##_base_ul((b), lep.up, had.ul); \
			lu = name##_base_lu((b), lep.lu, had.lu); \
			lp.y = name##_base_lt2((b), lep.lu, had.lt); \
			lp.z = name##_base_ll((b), lep.lp, had.ll); \
		} \
		break; \
	case 12: /* 1100 */ \
		{ \
			Lep##LP lep((kin)); \
			Had##LT had((kin), (sf)); \
			uu = name##_base_uu((b), lep.uu, had.uu); \
			up.x = name##_base_ut1((b), lep.up, had.ut); \
			lu = name##_base_lu((b), lep.lu, had.lu); \
			lp.x = name##_base_lt1((b), lep.lp, had.lt); \
		} \
		break; \
	case 13: /* 1101 */ \
		{ \
			Lep##LP lep((kin)); \
			Had##LP had((kin), (sf)); \
			uu = name##_base_uu((b), lep.uu, had.uu); \
			up.x = name##_base_ut1((b), lep.up, had.ut); \
			up.z = name##_base_ul((b), lep.up, had.ul); \
			lu = name##_base_lu((b), lep.lu, had.lu); \
			lp.x = name##_base_lt1((b), lep.lp, had.lt); \
			lp.z = name##_base_ll((b), lep.lp, had.ll); \
		} \
		break; \
	case 14: /* 1110 */ \
		{ \
			Lep##LP lep((kin)); \
			Had##LT had((kin), (sf)); \
			uu = name##_base_uu((b), lep.uu, had.uu); \
			up.x = name##_base_ut1((b), lep.up, had.ut); \
			up.y = name##_base_ut2((b), lep.uu, had.ut); \
			lu = name##_base_lu((b), lep.lu, had.lu); \
			lp.x = name##_base_lt1((b), lep.lp, had.lt); \
			lp.y = name##_base_lt2((b), lep.lu, had.lt); \
		} \
		break; \
	case 15: /* 1111 */ \
		{ \
			Lep##LP lep((kin)); \
			Had##LP had((kin), (sf)); \
			uu = name##_base_uu((b), lep.uu, had.uu); \
			up.x = name##_base_ut1((b), lep.up, had.ut); \
			up.y = name##_base_ut2((b), lep.uu, had.ut); \
			up.z = name##_base_ul((b), lep.up, had.ul); \
			lu = name##_base_lu((b), lep.lu, had.lu); \
			lp.x = name##_base_lt1((b), lep.lp, had.lt); \
			lp.y = name##_base_lt2((b), lep.lu, had.lt); \
			lp.z = name##_base_ll((b), lep.lp, had.ll); \
		} \
		break; \
	} \
	return uu + dot(up, (eta)) + (lambda_e)*(lu + dot(lp, (eta))); \
}())

// Similar to `SIDIS_MACRO_XS_FROM_BASE`, except this one works with base cross-
// sections where the XL, XT1, and XT2 cases are all grouped together into an
// XP case.
#define SIDIS_MACRO_XS_FROM_BASE_P(name, Lep, Had, kin, sf, b, lambda_e, eta) ([&]() { \
	/* Create a mask describing the polarization state. */ \
	unsigned pol_mask = (((lambda_e) != 0.) << 1) \
		| (((eta).x != 0. || (eta).y != 0. || (eta).z != 0.) << 0); \
	Real uu = 0.; \
	Vec3 up = VEC3_ZERO; \
	Real lu = 0.; \
	Vec3 lp = VEC3_ZERO; \
	switch (pol_mask) { \
	case 0:  /* 00 */ \
		{ \
			Lep##UU lep((kin)); \
			Had##UU had((kin), (sf)); \
			uu = name##_base_uu((b), lep.uu, had.uu); \
		} \
		break; \
	case 1:  /* 01 */ \
		{ \
			Lep##UP lep((kin)); \
			Had##UP had((kin), (sf)); \
			uu = name##_base_uu((b), lep.uu, had.uu); \
			up = name##_base_up((b), lep.uu, lep.up, had.up); \
		} \
		break; \
	case 2:  /* 10 */ \
		{ \
			Lep##LU lep((kin)); \
			Had##LU had((kin), (sf)); \
			uu = name##_base_uu((b), lep.uu, had.uu); \
			lu = name##_base_lu((b), lep.lu, had.lu); \
		} \
		break; \
	case 3:  /* 11 */ \
		{ \
			Lep##LP lep((kin)); \
			Had##LP had((kin), (sf)); \
			uu = name##_base_uu((b), lep.uu, had.uu); \
			up = name##_base_up((b), lep.uu, lep.up, had.up); \
			lu = name##_base_lu((b), lep.lu, had.lu); \
			lp = name##_base_lp((b), lep.lu, lep.lp, had.lp); \
		} \
		break; \
	} \
	return uu + dot(up, (eta)) + (lambda_e)*(lu + dot(lp, (eta))); \
}())

// This variant of `SIDIS_MACRO_XS_FROM_BASE_P` allows for an "endpoint" set of
// structure functions to be provided, for endpoint-subtraction-related
// calculations.
#define SIDIS_MACRO_XS_FROM_BASE_P_0(name, Lep, Had, kin, sf, had_0, b, lambda_e, eta) ([&]() { \
	/* Create a mask describing the polarization state. */ \
	unsigned pol_mask = (((lambda_e) != 0.) << 1) \
		| (((eta).x != 0. || (eta).y != 0. || (eta).z != 0.) << 0); \
	Real uu = 0.; \
	Vec3 up = VEC3_ZERO; \
	Real lu = 0.; \
	Vec3 lp = VEC3_ZERO; \
	switch (pol_mask) { \
	case 0:  /* 00 */ \
		{ \
			Lep##UU lep((kin)); \
			HadUU had_0_uu; \
			had_0_uu.uu = (had_0).uu; \
			Had##UU had((kin), (sf), had_0_uu); \
			uu = name##_base_uu((b), lep.uu, had.uu); \
		} \
		break; \
	case 1:  /* 01 */ \
		{ \
			Lep##UP lep((kin)); \
			HadUP had_0_up; \
			had_0_up.uu = (had_0).uu; \
			had_0_up.ul = (had_0).ul; \
			had_0_up.ut = (had_0).ut; \
			Had##UP had((kin), (sf), had_0_up); \
			uu = name##_base_uu((b), lep.uu, had.uu); \
			up = name##_base_up((b), lep.uu, lep.up, had.up); \
		} \
		break; \
	case 2:  /* 10 */ \
		{ \
			Lep##LU lep((kin)); \
			HadLU had_0_lu; \
			had_0_lu.uu = (had_0).uu; \
			had_0_lu.lu = (had_0).lu; \
			Had##LU had((kin), (sf), had_0_lu); \
			uu = name##_base_uu((b), lep.uu, had.uu); \
			lu = name##_base_lu((b), lep.lu, had.lu); \
		} \
		break; \
	case 3:  /* 11 */ \
		{ \
			Lep##LP lep((kin)); \
			Had##LP had((kin), (sf), (had_0)); \
			uu = name##_base_uu((b), lep.uu, had.uu); \
			up = name##_base_up((b), lep.uu, lep.up, had.up); \
			lu = name##_base_lu((b), lep.lu, had.lu); \
			lp = name##_base_lp((b), lep.lu, lep.lp, had.lp); \
		} \
		break; \
	} \
	return uu + dot(up, (eta)) + (lambda_e)*(lu + dot(lp, (eta))); \
}())

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

Real xs::born(Kinematics const& kin, Phenom const& phenom, SfSet const& sf, Real lambda_e, Vec3 eta) {
	Born b(kin, phenom);
	return SIDIS_MACRO_XS_FROM_BASE(born, LepBorn, Had, kin, sf, b, lambda_e, eta);
}

Real xs::amm(Kinematics const& kin, Phenom const& phenom, SfSet const& sf, Real lambda_e, Vec3 eta) {
	Amm b(kin, phenom);
	return SIDIS_MACRO_XS_FROM_BASE(amm, LepAmm, Had, kin, sf, b, lambda_e, eta);
}

Real xs::nrad_ir(Kinematics const& kin, Phenom const& phenom, SfSet const& sf, Real lambda_e, Vec3 eta, Real k_0_bar) {
	Nrad b(kin, phenom, k_0_bar);
	return SIDIS_MACRO_XS_FROM_BASE(nrad_ir, LepNrad, Had, kin, sf, b, lambda_e, eta);
}

Real xs::rad(KinematicsRad const& kin, Phenom const& phenom, SfSet const& sf, Real lambda_e, Vec3 eta) {
	Rad b(kin, phenom);
	return SIDIS_MACRO_XS_FROM_BASE_P(rad, LepRad, HadRad, kin, sf, b, lambda_e, eta);
}

Real xs::rad_f(KinematicsRad const& kin, Phenom const& phenom, SfSet const& sf, Real lambda_e, Vec3 eta) {
	Rad b(kin, phenom);
	return SIDIS_MACRO_XS_FROM_BASE_P(rad_f, LepRad, HadRadF, kin, sf, b, lambda_e, eta);
}

EstErr xs::nrad_integ(Kinematics const& kin, Phenom const& phenom, SfSet const& sf, Real lambda_e, Vec3 eta, Real k_0_bar, IntegParams params) {
	// The soft part of the radiative cross-section (below `k_0_bar`) is bundled
	// into the return value here.
	Real xs_nrad_ir = nrad_ir(kin, phenom, sf, lambda_e, eta, k_0_bar);
	// TODO: The integration parameters should be modified here to account for
	// the `xs_nrad_ir` contribution.
	EstErr xs_rad_f = rad_f_integ(kin, phenom, sf, lambda_e, eta, k_0_bar, params);
	return { xs_nrad_ir + xs_rad_f.val, xs_rad_f.err };
}

EstErr xs::rad_f_integ(Kinematics const& kin, Phenom const& phenom, SfSet const& sf, Real lambda_e, Vec3 eta, Real k_0_bar, IntegParams params) {
	HadLP had_0(kin, sf);
	CutRad cut;
	cut.k_0_bar = Bound(0., k_0_bar);
	EstErr xs_integ = integrate<3>(
		[&](std::array<Real, 3> x) {
			KinematicsRad kin_rad;
			Real jac;
			if (!take(cut, kin, x.data(), &kin_rad, &jac)) {
				return 0.;
			}
			Rad b(kin_rad, phenom);
			Real xs = jac * SIDIS_MACRO_XS_FROM_BASE_P_0(rad_f, LepRad, HadRadF, kin_rad, sf, had_0, b, lambda_e, eta);
			if (std::isnan(xs)) {
				// If the result is `NaN`, it most likely means we went out of
				// the allowed region for the structure function grids (or we
				// are in a kinematically disallowed region). In that case, just
				// return zero.
				return 0.;
			} else {
				return xs;
			}
		},
		std::array<Real, 3>{ 0., 0., 0. },
		std::array<Real, 3>{ 1., 1., 1. },
		params);
	return xs_integ;
}

EstErr xs::rad_integ(Kinematics const& kin, Phenom const& phenom, SfSet const& sf, Real lambda_e, Vec3 eta, Real k_0_bar, IntegParams params) {
	CutRad cut;
	cut.k_0_bar = Bound(k_0_bar, INF);
	EstErr xs_integ = integrate<3>(
		[&](std::array<Real, 3> x) {
			KinematicsRad kin_rad;
			Real jac;
			if (!take(cut, kin, x.data(), &kin_rad, &jac)) {
				return 0.;
			}
			Rad b(kin_rad, phenom);
			Real xs = jac * SIDIS_MACRO_XS_FROM_BASE_P(rad, LepRad, HadRad, kin_rad, sf, b, lambda_e, eta);
			if (std::isnan(xs)) {
				return 0.;
			} else {
				return xs;
			}
		},
		std::array<Real, 3>{ 0., 0., 0. },
		std::array<Real, 3>{ 1., 1., 1. },
		params);
	return xs_integ;
}

Real xs::born(Kinematics const& kin, SfSet const& sf, Real lambda_e, Vec3 eta) {
	return born(kin, Phenom(kin), sf, lambda_e, eta);
}
Real xs::amm(Kinematics const& kin, SfSet const& sf, Real lambda_e, Vec3 eta) {
	return amm(kin, Phenom(kin), sf, lambda_e, eta);
}
Real xs::nrad_ir(Kinematics const& kin, SfSet const& sf, Real lambda_e, Vec3 eta, Real k_0_bar) {
	return nrad_ir(kin, Phenom(kin), sf, lambda_e, eta, k_0_bar);
}
Real xs::rad(KinematicsRad const& kin, SfSet const& sf, Real lambda_e, Vec3 eta) {
	return rad(kin, Phenom(kin.project()), sf, lambda_e, eta);
}
Real xs::rad_f(KinematicsRad const& kin, SfSet const& sf, Real lambda_e, Vec3 eta) {
	return rad_f(kin, Phenom(kin.project()), sf, lambda_e, eta);
}
EstErr xs::nrad_integ(Kinematics const& kin, SfSet const& sf, Real lambda_e, Vec3 eta, Real k_0_bar, IntegParams params) {
	return nrad_integ(kin, Phenom(kin), sf, lambda_e, eta, k_0_bar, params);
}
EstErr xs::rad_f_integ(Kinematics const& kin, SfSet const& sf, Real lambda_e, Vec3 eta, Real k_0_bar, IntegParams params) {
	return rad_f_integ(kin, Phenom(kin), sf, lambda_e, eta, k_0_bar, params);
}
EstErr xs::rad_integ(Kinematics const& kin, SfSet const& sf, Real lambda_e, Vec3 eta, Real k_0_bar, IntegParams params) {
	return rad_integ(kin, Phenom(kin), sf, lambda_e, eta, k_0_bar, params);
}

// Radiative corrections to Born cross-section.
Real xs::delta_vert_rad_ir(Kinematics const& kin, Real k_0_bar) {
	// Paragraph following equation [1.C17].
	Real k0_max = (kin.mx_sq - sq(kin.Mth))/(2.*kin.mx);
	if (!(k_0_bar > 0.)) {
		return -INF;
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
		return INF;
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

// Born base functions.
Born::Born(Kinematics const& kin, Phenom const& phenom) :
	// Equation [1.15]. The `Q^4` factor has been absorbed into `C_1`.
	coeff((sq(phenom.alpha_qed)*kin.S*sq(kin.S_x))/(8.*kin.M*kin.ph_l*kin.lambda_S)) { }

Real xs::born_base_uu(Born const& b, LepBornBaseUU const& lep, HadBaseUU const& had) {
	return b.coeff*(
		lep.theta_1*had.H_10
		+ lep.theta_2*had.H_20
		+ lep.theta_3*had.H_30
		+ lep.theta_4*had.H_40);
}
Real xs::born_base_ul(Born const& b, LepBornBaseUP const& lep, HadBaseUL const& had) {
	return b.coeff*(lep.theta_6*had.H_63 + lep.theta_8*had.H_83);
}
Real xs::born_base_ut1(Born const& b, LepBornBaseUP const& lep, HadBaseUT const& had) {
	return b.coeff*(lep.theta_6*had.H_61 + lep.theta_8*had.H_81);
}
Real xs::born_base_ut2(Born const& b, LepBornBaseUU const& lep, HadBaseUT const& had) {
	return b.coeff*(
		lep.theta_1*had.H_12
		+ lep.theta_2*had.H_22
		+ lep.theta_3*had.H_32
		+ lep.theta_4*had.H_42);
}
Real xs::born_base_lu(Born const& b, LepBornBaseLU const& lep, HadBaseLU const& had) {
	return b.coeff*lep.theta_5*had.H_50;
}
Real xs::born_base_ll(Born const& b, LepBornBaseLP const& lep, HadBaseLL const& had) {
	return b.coeff*(lep.theta_7*had.H_73 + lep.theta_9*had.H_93);
}
Real xs::born_base_lt1(Born const& b, LepBornBaseLP const& lep, HadBaseLT const& had) {
	return b.coeff*(lep.theta_7*had.H_71 + lep.theta_9*had.H_91);
}
Real xs::born_base_lt2(Born const& b, LepBornBaseLU const& lep, HadBaseLT const& had) {
	return b.coeff*lep.theta_5*had.H_52;
}

// AMM base functions.
Amm::Amm(Kinematics const& kin, Phenom const& phenom) {
	// Equation [1.53]. The `Q^4` factor has been absorbed into `C_1`.
	Real lambda_m = kin.Q_sq*(kin.Q_sq + 4.*sq(kin.m));
	Real lambda_m_sqrt = std::sqrt(lambda_m);
	Real diff_m = sqrt1p_1m((4.*sq(kin.m))/kin.Q_sq);
	Real sum_m = 2. + diff_m;
	Real L_m = 1./lambda_m_sqrt*std::log(sum_m/diff_m);
	coeff = L_m*kin.Q_sq*(std::pow(phenom.alpha_qed, 3)*sq(kin.m)*kin.S*sq(kin.S_x))
		/(16.*PI*kin.M*kin.ph_l*kin.lambda_S);
}

Real xs::amm_base_uu(Amm const& b, LepAmmBaseUU const& lep, HadBaseUU const& had) {
	return b.coeff*(
		lep.theta_1*had.H_10
		+ lep.theta_2*had.H_20
		+ lep.theta_3*had.H_30
		+ lep.theta_4*had.H_40);
}
Real xs::amm_base_ul(Amm const& b, LepAmmBaseUP const& lep, HadBaseUL const& had) {
	return b.coeff*(lep.theta_6*had.H_63 + lep.theta_8*had.H_83);
}
Real xs::amm_base_ut1(Amm const& b, LepAmmBaseUP const& lep, HadBaseUT const& had) {
	return b.coeff*(lep.theta_6*had.H_61 + lep.theta_8*had.H_81);
}
Real xs::amm_base_ut2(Amm const& b, LepAmmBaseUU const& lep, HadBaseUT const& had) {
	return b.coeff*(
		lep.theta_1*had.H_12
		+ lep.theta_2*had.H_22
		+ lep.theta_3*had.H_32
		+ lep.theta_4*had.H_42);
}
Real xs::amm_base_lu(Amm const& b, LepAmmBaseLU const& lep, HadBaseLU const& had) {
	return b.coeff*lep.theta_5*had.H_50;
}
Real xs::amm_base_ll(Amm const& b, LepAmmBaseLP const& lep, HadBaseLL const& had) {
	return b.coeff*(lep.theta_7*had.H_73 + lep.theta_9*had.H_93);
}
Real xs::amm_base_lt1(Amm const& b, LepAmmBaseLP const& lep, HadBaseLT const& had) {
	return b.coeff*(lep.theta_7*had.H_71 + lep.theta_9*had.H_91);
}
Real xs::amm_base_lt2(Amm const& b, LepAmmBaseLU const& lep, HadBaseLT const& had) {
	return b.coeff*lep.theta_5*had.H_52;
}

// Non-radiative infrared-divergence-free base functions.
Nrad::Nrad(Kinematics const& kin, Phenom const& phenom, Real k_0_bar) {
	Born born(kin, phenom);
	Amm amm(kin, phenom);
	Real born_factor = 1. + phenom.alpha_qed/PI*(
		delta_vert_rad_ir(kin, k_0_bar)
		+ delta_vac_lep(kin)
		+ phenom.delta_vac_had);
	coeff_born = born_factor*born.coeff;
	coeff_amm = amm.coeff;
}

Real xs::nrad_ir_base_uu(Nrad const& b, LepNradBaseUU const& lep, HadBaseUU const& had) {
	return
		(b.coeff_born*lep.born.theta_1 + b.coeff_amm*lep.amm.theta_1)*had.H_10
		+ (b.coeff_born*lep.born.theta_2 + b.coeff_amm*lep.amm.theta_2)*had.H_20
		+ (b.coeff_born*lep.born.theta_3 + b.coeff_amm*lep.amm.theta_3)*had.H_30
		+ (b.coeff_born*lep.born.theta_4 + b.coeff_amm*lep.amm.theta_4)*had.H_40;
}
Real xs::nrad_ir_base_ul(Nrad const& b, LepNradBaseUP const& lep, HadBaseUL const& had) {
	return
		(b.coeff_born*lep.born.theta_6 + b.coeff_amm*lep.amm.theta_6)*had.H_63
		+ (b.coeff_born*lep.born.theta_8 + b.coeff_amm*lep.amm.theta_8)*had.H_83;
}
Real xs::nrad_ir_base_ut1(Nrad const& b, LepNradBaseUP const& lep, HadBaseUT const& had) {
	return
		(b.coeff_born*lep.born.theta_6 + b.coeff_amm*lep.amm.theta_6)*had.H_61
		+ (b.coeff_born*lep.born.theta_8 + b.coeff_amm*lep.amm.theta_8)*had.H_81;
}
Real xs::nrad_ir_base_ut2(Nrad const& b, LepNradBaseUU const& lep, HadBaseUT const& had) {
	return
		(b.coeff_born*lep.born.theta_1 + b.coeff_amm*lep.amm.theta_1)*had.H_12
		+ (b.coeff_born*lep.born.theta_2 + b.coeff_amm*lep.amm.theta_2)*had.H_22
		+ (b.coeff_born*lep.born.theta_3 + b.coeff_amm*lep.amm.theta_3)*had.H_32
		+ (b.coeff_born*lep.born.theta_4 + b.coeff_amm*lep.amm.theta_4)*had.H_42;
}
Real xs::nrad_ir_base_lu(Nrad const& b, LepNradBaseLU const& lep, HadBaseLU const& had) {
	return (b.coeff_born*lep.born.theta_5 + b.coeff_amm*lep.amm.theta_5)*had.H_50;
}
Real xs::nrad_ir_base_ll(Nrad const& b, LepNradBaseLP const& lep, HadBaseLL const& had) {
	return
		(b.coeff_born*lep.born.theta_7 + b.coeff_amm*lep.amm.theta_7)*had.H_73
		+ (b.coeff_born*lep.born.theta_9 + b.coeff_amm*lep.amm.theta_9)*had.H_93;
}
Real xs::nrad_ir_base_lt1(Nrad const& b, LepNradBaseLP const& lep, HadBaseLT const& had) {
	return
		(b.coeff_born*lep.born.theta_7 + b.coeff_amm*lep.amm.theta_7)*had.H_71
		+ (b.coeff_born*lep.born.theta_9 + b.coeff_amm*lep.amm.theta_9)*had.H_91;
}
Real xs::nrad_ir_base_lt2(Nrad const& b, LepNradBaseLU const& lep, HadBaseLT const& had) {
	return (b.coeff_born*lep.born.theta_5 + b.coeff_amm*lep.amm.theta_5)*had.H_52;
}

// Radiative base functions.
Rad::Rad(KinematicsRad const& kin, Phenom const& phenom) {
	// Equation [1.43].
	coeff = -(std::pow(phenom.alpha_qed, 3)*kin.S*sq(kin.S_x))
		/(64.*sq(PI)*kin.M*kin.ph_l*kin.lambda_S*kin.lambda_Y_sqrt);
	R = kin.R;
}

Real xs::rad_base_uu(Rad const& b, LepRadBaseUU const& lep, HadRadBaseUU const& had) {
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
Vec3 xs::rad_base_up(Rad const& b, LepRadBaseUU const& lep_uu, LepRadBaseUP const& lep_up, HadRadBaseUP const& had) {
	return b.coeff*(
		1./b.R*(
			lep_uu.theta_011*had.H_1
			+ lep_uu.theta_021*had.H_2
			+ lep_uu.theta_031*had.H_3
			+ lep_uu.theta_041*had.H_4
			+ lep_up.theta_061*had.H_6
			+ lep_up.theta_081*had.H_8)
		+ (
			lep_uu.theta_012*had.H_1
			+ lep_uu.theta_022*had.H_2
			+ lep_uu.theta_032*had.H_3
			+ lep_uu.theta_042*had.H_4
			+ lep_up.theta_062*had.H_6
			+ lep_up.theta_082*had.H_8)
		+ b.R*(
			lep_uu.theta_013*had.H_1
			+ lep_uu.theta_023*had.H_2
			+ lep_uu.theta_033*had.H_3
			+ lep_uu.theta_043*had.H_4
			+ lep_up.theta_063*had.H_6
			+ lep_up.theta_083*had.H_8)
		+ b.R*b.R*(
			lep_up.theta_064*had.H_6
			+ lep_up.theta_084*had.H_8));
}
Real xs::rad_base_lu(Rad const& b, LepRadBaseLU const& lep, HadRadBaseLU const& had) {
	return b.coeff*(
		1./b.R*(lep.theta_051 + lep.theta_151)*had.H_50
		+ (lep.theta_052 + lep.theta_152)*had.H_50
		+ b.R*(lep.theta_053 + lep.theta_153)*had.H_50);
}
Vec3 xs::rad_base_lp(Rad const& b, LepRadBaseLU const& lep_lu, LepRadBaseLP const& lep_lp, HadRadBaseLP const& had) {
	return b.coeff*(
		1./b.R*(
			(lep_lu.theta_051 + lep_lu.theta_151)*had.H_5
			+ (lep_lp.theta_071 + lep_lp.theta_171)*had.H_7
			+ (lep_lp.theta_091 + lep_lp.theta_191)*had.H_9)
		+ (
		 	(lep_lu.theta_052 + lep_lu.theta_152)*had.H_5
			+ (lep_lp.theta_072 + lep_lp.theta_172)*had.H_7
			+ (lep_lp.theta_092 + lep_lp.theta_192)*had.H_9)
		+ b.R*(
		 	(lep_lu.theta_053 + lep_lu.theta_153)*had.H_5
			+ (lep_lp.theta_073 + lep_lp.theta_173)*had.H_7
			+ (lep_lp.theta_093 + lep_lp.theta_193)*had.H_9)
		+ b.R*b.R*(
			(lep_lp.theta_074 + lep_lp.theta_174)*had.H_7
			+ (lep_lp.theta_094 + lep_lp.theta_194)*had.H_9));
}

Real xs::rad_f_base_uu(Rad const& b, LepRadBaseUU const& lep, HadRadFBaseUU const& had) {
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
Vec3 xs::rad_f_base_up(Rad const& b, LepRadBaseUU const& lep_uu, LepRadBaseUP const& lep_up, HadRadFBaseUP const& had) {
	return b.coeff*(
		(
			lep_uu.theta_011*had.H_1_diff
			+ lep_uu.theta_021*had.H_2_diff
			+ lep_uu.theta_031*had.H_3_diff
			+ lep_uu.theta_041*had.H_4_diff
			+ lep_up.theta_061*had.H_6_diff
			+ lep_up.theta_081*had.H_8_diff)
		+ (
			lep_uu.theta_012*had.H_1
			+ lep_uu.theta_022*had.H_2
			+ lep_uu.theta_032*had.H_3
			+ lep_uu.theta_042*had.H_4
			+ lep_up.theta_062*had.H_6
			+ lep_up.theta_082*had.H_8)
		+ b.R*(
			lep_uu.theta_013*had.H_1
			+ lep_uu.theta_023*had.H_2
			+ lep_uu.theta_033*had.H_3
			+ lep_uu.theta_043*had.H_4
			+ lep_up.theta_063*had.H_6
			+ lep_up.theta_083*had.H_8)
		+ b.R*b.R*(
			lep_up.theta_064*had.H_6
			+ lep_up.theta_084*had.H_8));
}
Real xs::rad_f_base_lu(Rad const& b, LepRadBaseLU const& lep, HadRadFBaseLU const& had) {
	return b.coeff*(
		(lep.theta_051 + lep.theta_151)*had.H_50_diff
		+ (lep.theta_052 + lep.theta_152)*had.H_50
		+ b.R*(lep.theta_053 + lep.theta_153)*had.H_50);
}
Vec3 xs::rad_f_base_lp(Rad const& b, LepRadBaseLU const& lep_lu, LepRadBaseLP const& lep_lp, HadRadFBaseLP const& had) {
	return b.coeff*(
		(
			(lep_lu.theta_051 + lep_lu.theta_151)*had.H_5_diff
			+ (lep_lp.theta_071 + lep_lp.theta_171)*had.H_7_diff
			+ (lep_lp.theta_091 + lep_lp.theta_191)*had.H_9_diff)
		+ (
		 	(lep_lu.theta_052 + lep_lu.theta_152)*had.H_5
			+ (lep_lp.theta_072 + lep_lp.theta_172)*had.H_7
			+ (lep_lp.theta_092 + lep_lp.theta_192)*had.H_9)
		+ b.R*(
		 	(lep_lu.theta_053 + lep_lu.theta_153)*had.H_5
			+ (lep_lp.theta_073 + lep_lp.theta_173)*had.H_7
			+ (lep_lp.theta_093 + lep_lp.theta_193)*had.H_9)
		+ b.R*b.R*(
			(lep_lp.theta_074 + lep_lp.theta_174)*had.H_7
			+ (lep_lp.theta_094 + lep_lp.theta_194)*had.H_9));
}

