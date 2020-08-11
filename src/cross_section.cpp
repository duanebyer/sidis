#include "sidis/cross_section.hpp"

#include <cmath>

#include "sidis/constant.hpp"
#include "sidis/hadronic_coeff.hpp"
#include "sidis/kinematics.hpp"
#include "sidis/leptonic_coeff.hpp"
#include "sidis/math.hpp"
#include "sidis/structure_function.hpp"
#include "sidis/vector.hpp"

using namespace sidis;
using namespace sidis::xs;
using namespace sidis::constant;
using namespace sidis::had;
using namespace sidis::kin;
using namespace sidis::lep;
using namespace sidis::math;
using namespace sidis::sf;

namespace {

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
	// Equation [1.15].
	coeff = (sq(ALPHA)*kin.S*sq(kin.S_x))
		/(8.*kin.M*sq(kin.Q_sq)*kin.ph_l*kin.lambda_S);
	// Equation [1.18].
	c_1 = (4.*kin.M*kin.ph_l*(kin.Q_sq + 2.*kin.x*sq(kin.M)))/sq(kin.Q_sq);
}

Amm::Amm(Kinematics kin) {
	// Equation [1.53].
	coeff = (std::pow(ALPHA, 3)*sq(kin.m)*kin.S*sq(kin.S_x))
		/(16.*PI*kin.M*kin.Q_sq*kin.ph_l*kin.lambda_S);
	// Equation [1.18].
	c_1 = (4.*kin.M*kin.ph_l*(kin.Q_sq + 2.*kin.x*sq(kin.M)))/sq(kin.Q_sq);
}

Real xs::born(Real lambda_e, Vec3 eta, Kinematics kin, Sf sf) {
	// Calculate the complete Born cross-section by interpolating all of the
	// individual pieces together.
	Born b(kin);

	LepBornUU lep_uu(kin);
	LepBornUP lep_up(kin);
	HadUU had_uu(kin, sf);
	HadUL had_ul(kin, sf);
	HadUT1 had_ut1(kin, sf);
	HadUT2 had_ut2(kin, sf);

	LepBornLU lep_lu(kin);
	LepBornLP lep_lp(kin);
	HadLU had_lu(kin, sf);
	HadLL had_ll(kin, sf);
	HadLT1 had_lt1(kin, sf);
	HadLT2 had_lt2(kin, sf);

	Real uu = born_uu(b, lep_uu, had_uu);
	Real ul = born_ul(b, lep_up, had_ul);
	Real ut1 = born_ut1(b, lep_up, had_ut1);
	Real ut2 = born_ut2(b, lep_uu, had_ut2);
	Real lu = born_lu(b, lep_lu, had_lu);
	Real ll = born_ll(b, lep_lp, had_ll);
	Real lt1 = born_lt1(b, lep_lp, had_lt1);
	Real lt2 = born_lt2(b, lep_lu, had_lt2);

	Real up = dot(eta, Vec3(ut1, ut2, ul));
	Real lp = dot(eta, Vec3(lt1, lt2, ll));
	return uu + up + lambda_e * (lu + lp);
}

Real xs::amm(Real lambda_e, Vec3 eta, Kinematics kin, Sf sf) {
	Amm b(kin);

	LepAmmUU lep_uu(kin);
	LepAmmUP lep_up(kin);
	HadUU had_uu(kin, sf);
	HadUL had_ul(kin, sf);
	HadUT1 had_ut1(kin, sf);
	HadUT2 had_ut2(kin, sf);

	LepAmmLU lep_lu(kin);
	LepAmmLP lep_lp(kin);
	HadLU had_lu(kin, sf);
	HadLL had_ll(kin, sf);
	HadLT1 had_lt1(kin, sf);
	HadLT2 had_lt2(kin, sf);

	Real uu = amm_uu(b, lep_uu, had_uu);
	Real ul = amm_ul(b, lep_up, had_ul);
	Real ut1 = amm_ut1(b, lep_up, had_ut1);
	Real ut2 = amm_ut2(b, lep_uu, had_ut2);
	Real lu = amm_lu(b, lep_lu, had_lu);
	Real ll = amm_ll(b, lep_lp, had_ll);
	Real lt1 = amm_lt1(b, lep_lp, had_lt1);
	Real lt2 = amm_lt2(b, lep_lu, had_lt2);

	Real up = dot(eta, Vec3(ut1, ut2, ul));
	Real lp = dot(eta, Vec3(lt1, lt2, ll));
	return uu + up + lambda_e * (lu + lp);
}

// Born base functions.
Real xs::born_uu(Born b, LepBornUU lep, HadUU had) {
	return b.coeff*b.c_1*(
		lep.theta_1*had.H_1
		+ lep.theta_2*had.H_2
		+ lep.theta_3*had.H_3
		+ lep.theta_4*had.H_4);
}
Real xs::born_ul(Born b, LepBornUP lep, HadUL had) {
	return b.coeff*b.c_1*(lep.theta_6*had.H_6 + lep.theta_8*had.H_8);
}
Real xs::born_ut1(Born b, LepBornUP lep, HadUT1 had) {
	return b.coeff*b.c_1*(lep.theta_6*had.H_6 + lep.theta_8*had.H_8);
}
Real xs::born_ut2(Born b, LepBornUU lep, HadUT2 had) {
	return b.coeff*b.c_1*(
		lep.theta_1*had.H_1
		+ lep.theta_2*had.H_2
		+ lep.theta_3*had.H_3
		+ lep.theta_4*had.H_4);
}
Real xs::born_lu(Born b, LepBornLU lep, HadLU had) {
	return b.coeff*b.c_1*lep.theta_5*had.H_5;
}
Real xs::born_ll(Born b, LepBornLP lep, HadLL had) {
	return b.coeff*b.c_1*(lep.theta_7*had.H_7 + lep.theta_9*had.H_9);
}
Real xs::born_lt1(Born b, LepBornLP lep, HadLT1 had) {
	return b.coeff*b.c_1*(lep.theta_7*had.H_7 + lep.theta_9*had.H_9);
}
Real xs::born_lt2(Born b, LepBornLU lep, HadLT2 had) {
	return b.coeff*b.c_1*lep.theta_5*had.H_5;
}

// AMM base functions.
Real xs::amm_uu(Amm b, LepAmmUU lep, HadUU had) {
	return b.coeff*b.c_1*(
		lep.theta_1*had.H_1
		+ lep.theta_2*had.H_2
		+ lep.theta_3*had.H_3
		+ lep.theta_4*had.H_4);
}
Real xs::amm_ul(Amm b, LepAmmUP lep, HadUL had) {
	return b.coeff*b.c_1*(lep.theta_6*had.H_6 + lep.theta_8*had.H_8);
}
Real xs::amm_ut1(Amm b, LepAmmUP lep, HadUT1 had) {
	return b.coeff*b.c_1*(lep.theta_6*had.H_6 + lep.theta_8*had.H_8);
}
Real xs::amm_ut2(Amm b, LepAmmUU lep, HadUT2 had) {
	return b.coeff*b.c_1*(
		lep.theta_1*had.H_1
		+ lep.theta_2*had.H_2
		+ lep.theta_3*had.H_3
		+ lep.theta_4*had.H_4);
}
Real xs::amm_lu(Amm b, LepAmmLU lep, HadLU had) {
	return b.coeff*b.c_1*lep.theta_5*had.H_5;
}
Real xs::amm_ll(Amm b, LepAmmLP lep, HadLL had) {
	return b.coeff*b.c_1*(lep.theta_7*had.H_7 + lep.theta_9*had.H_9);
}
Real xs::amm_lt1(Amm b, LepAmmLP lep, HadLT1 had) {
	return b.coeff*b.c_1*(lep.theta_7*had.H_7 + lep.theta_9*had.H_9);
}
Real xs::amm_lt2(Amm b, LepAmmLU lep, HadLT2 had) {
	return b.coeff*b.c_1*lep.theta_5*had.H_5;
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
	Real d = lambda_X_prime*diff_m;
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
	Real z_1d = 1./lambda_X_prime_sqrt*(
		X_prime*sum_X_prime*(S_prime - X_prime) - c_diff*(d + d_diff));
	Real z_2d = 1./lambda_X_prime_sqrt*(
		X_prime*sum_X_prime*(S_prime - X_prime) + c_diff*(d + d_sum));
	Real z_3d = 1./lambda_X_prime_sqrt*(
		-X_prime*diff_X_prime*(S_prime - X_prime) + c_sum*(d + d_diff));
	Real z_4d = 1./lambda_X_prime_sqrt*(
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
	Real S_phi_a = -Q_m_sq/(2.*lambda_m_sqrt);
	Real S_phi_b = std::log(
		(X_prime - lambda_X_prime_sqrt)/(X_prime + lambda_X_prime_sqrt));
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

Real xs::delta_vac_lep(kin::Kinematics kin) {
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

Real xs::delta_vac_had(kin::Kinematics kin) {
	if (kin.Q_sq < 1.) {
		return -(2.*PI)/ALPHA*(-1.345e-9 - 2.302e-3*std::log(1. + 4.091*kin.Q_sq));
	} else if (kin.Q_sq < 64.) {
		return -(2.*PI)/ALPHA*(-1.512e-3 - 2.822e-3*std::log(1. + 1.218*kin.Q_sq));
	} else {
		return -(2.*PI)/ALPHA*(-1.1344-3 - 3.0680-3*std::log(1. + 0.99992*kin.Q_sq));
	}
}

