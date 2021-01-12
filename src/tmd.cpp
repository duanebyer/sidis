#include "sidis/tmd.hpp"

#include <cmath>

#include "sidis/constant.hpp"

using namespace sidis;
using namespace sidis::constant;
using namespace sidis::sf;

namespace {

Real gaussian(Real k_perp_sq, Real mean) {
	return std::exp(-k_perp_sq/mean)/(PI*mean);
}

}

// Standard TMDs.
Real TmdSet::xf1(unsigned, Real, Real, Real) const {
	return 0.;
}
Real TmdSet::xf1Tperp(unsigned, Real, Real, Real) const {
	return 0.;
}
Real TmdSet::xfT(unsigned, Real, Real, Real) const {
	return 0.;
}
Real TmdSet::xfperp(unsigned, Real, Real, Real) const {
	return 0.;
}
Real TmdSet::xfLperp(unsigned, Real, Real, Real) const {
	return 0.;
}
Real TmdSet::xfTperp(unsigned, Real, Real, Real) const {
	return 0.;
}
Real TmdSet::xg1(unsigned, Real, Real, Real) const {
	return 0.;
}
Real TmdSet::xg1Tperp(unsigned, Real, Real, Real) const {
	return 0.;
}
Real TmdSet::xgT(unsigned, Real, Real, Real) const {
	return 0.;
}
Real TmdSet::xgperp(unsigned, Real, Real, Real) const {
	return 0.;
}
Real TmdSet::xgLperp(unsigned, Real, Real, Real) const {
	return 0.;
}
Real TmdSet::xgTperp(unsigned, Real, Real, Real) const {
	return 0.;
}
Real TmdSet::xh1(unsigned, Real, Real, Real) const {
	return 0.;
}
Real TmdSet::xh1perp(unsigned, Real, Real, Real) const {
	return 0.;
}
Real TmdSet::xh1Lperp(unsigned, Real, Real, Real) const {
	return 0.;
}
Real TmdSet::xh1Tperp(unsigned, Real, Real, Real) const {
	return 0.;
}
Real TmdSet::xh(unsigned, Real, Real, Real) const {
	return 0.;
}
Real TmdSet::xhL(unsigned, Real, Real, Real) const {
	return 0.;
}
Real TmdSet::xhT(unsigned, Real, Real, Real) const {
	return 0.;
}
Real TmdSet::xhTperp(unsigned, Real, Real, Real) const {
	return 0.;
}
Real TmdSet::xe(unsigned, Real, Real, Real) const {
	return 0.;
}
Real TmdSet::xeL(unsigned, Real, Real, Real) const {
	return 0.;
}
Real TmdSet::xeT(unsigned, Real, Real, Real) const {
	return 0.;
}
Real TmdSet::xeTperp(unsigned, Real, Real, Real) const {
	return 0.;
}

Real TmdSet::D1(Hadron, unsigned, Real, Real, Real) const {
	return 0.;
}
Real TmdSet::H1perp(Hadron, unsigned, Real, Real, Real) const {
	return 0.;
}
Real TmdSet::Dperp_tilde(Hadron, unsigned, Real, Real, Real) const {
	return 0.;
}
Real TmdSet::H_tilde(Hadron, unsigned, Real, Real, Real) const {
	return 0.;
}
Real TmdSet::Gperp_tilde(Hadron, unsigned, Real, Real, Real) const {
	return 0.;
}
Real TmdSet::E_tilde(Hadron, unsigned, Real, Real, Real) const {
	return 0.;
}

// Gaussian approximation.
TmdGaussianSet::TmdGaussianSet(
	unsigned flavor_count,
	constant::Nucleus target,
	Real mean_f1,
	Real mean_f1Tperp,
	Real mean_fT,
	Real mean_fperp,
	Real mean_fLperp,
	Real mean_fTperp,
	Real mean_g1,
	Real mean_g1Tperp,
	Real mean_gT,
	Real mean_gperp,
	Real mean_gLperp,
	Real mean_gTperp,
	Real mean_h1,
	Real mean_h1perp,
	Real mean_h1Lperp,
	Real mean_h1Tperp,
	Real mean_h,
	Real mean_hL,
	Real mean_hT,
	Real mean_hTperp,
	Real mean_e,
	Real mean_eL,
	Real mean_eT,
	Real mean_eTperp,
	Real mean_D1,
	Real mean_H1perp,
	Real mean_Dperp_tilde,
	Real mean_H_tilde,
	Real mean_Gperp_tilde,
	Real mean_E_tilde) :
	TmdSet(flavor_count, target),
	mean_f1(mean_f1),
	mean_f1Tperp(mean_f1Tperp),
	mean_fT(mean_fT),
	mean_fperp(mean_fperp),
	mean_fLperp(mean_fLperp),
	mean_fTperp(mean_fTperp),
	mean_g1(mean_g1),
	mean_g1Tperp(mean_g1Tperp),
	mean_gT(mean_gT),
	mean_gperp(mean_gperp),
	mean_gLperp(mean_gLperp),
	mean_gTperp(mean_gTperp),
	mean_h1(mean_h1),
	mean_h1perp(mean_h1perp),
	mean_h1Lperp(mean_h1Lperp),
	mean_h1Tperp(mean_h1Tperp),
	mean_h(mean_h),
	mean_hL(mean_hL),
	mean_hT(mean_hT),
	mean_hTperp(mean_hTperp),
	mean_e(mean_e),
	mean_eL(mean_eL),
	mean_eT(mean_eT),
	mean_eTperp(mean_eTperp),
	mean_D1(mean_D1),
	mean_H1perp(mean_H1perp),
	mean_Dperp_tilde(mean_Dperp_tilde),
	mean_H_tilde(mean_H_tilde),
	mean_Gperp_tilde(mean_Gperp_tilde),
	mean_E_tilde(mean_E_tilde) { }

Real TmdGaussianSet::xf1(unsigned, Real, Real) const {
	return 0.;
}
Real TmdGaussianSet::xf1Tperp(unsigned, Real, Real) const {
	return 0.;
}
Real TmdGaussianSet::xfT(unsigned, Real, Real) const {
	return 0.;
}
Real TmdGaussianSet::xfperp(unsigned, Real, Real) const {
	return 0.;
}
Real TmdGaussianSet::xfLperp(unsigned, Real, Real) const {
	return 0.;
}
Real TmdGaussianSet::xfTperp(unsigned, Real, Real) const {
	return 0.;
}
Real TmdGaussianSet::xg1(unsigned, Real, Real) const {
	return 0.;
}
Real TmdGaussianSet::xg1Tperp(unsigned, Real, Real) const {
	return 0.;
}
Real TmdGaussianSet::xgT(unsigned, Real, Real) const {
	return 0.;
}
Real TmdGaussianSet::xgperp(unsigned, Real, Real) const {
	return 0.;
}
Real TmdGaussianSet::xgLperp(unsigned, Real, Real) const {
	return 0.;
}
Real TmdGaussianSet::xgTperp(unsigned, Real, Real) const {
	return 0.;
}
Real TmdGaussianSet::xh1(unsigned, Real, Real) const {
	return 0.;
}
Real TmdGaussianSet::xh1perp(unsigned, Real, Real) const {
	return 0.;
}
Real TmdGaussianSet::xh1Lperp(unsigned, Real, Real) const {
	return 0.;
}
Real TmdGaussianSet::xh1Tperp(unsigned, Real, Real) const {
	return 0.;
}
Real TmdGaussianSet::xh(unsigned, Real, Real) const {
	return 0.;
}
Real TmdGaussianSet::xhL(unsigned, Real, Real) const {
	return 0.;
}
Real TmdGaussianSet::xhT(unsigned, Real, Real) const {
	return 0.;
}
Real TmdGaussianSet::xhTperp(unsigned, Real, Real) const {
	return 0.;
}
Real TmdGaussianSet::xe(unsigned, Real, Real) const {
	return 0.;
}
Real TmdGaussianSet::xeL(unsigned, Real, Real) const {
	return 0.;
}
Real TmdGaussianSet::xeT(unsigned, Real, Real) const {
	return 0.;
}
Real TmdGaussianSet::xeTperp(unsigned, Real, Real) const {
	return 0.;
}

Real TmdGaussianSet::D1(Hadron, unsigned, Real, Real) const {
	return 0.;
}
Real TmdGaussianSet::H1perp(Hadron, unsigned, Real, Real) const {
	return 0.;
}
Real TmdGaussianSet::Dperp_tilde(Hadron, unsigned, Real, Real) const {
	return 0.;
}
Real TmdGaussianSet::H_tilde(Hadron, unsigned, Real, Real) const {
	return 0.;
}
Real TmdGaussianSet::Gperp_tilde(Hadron, unsigned, Real, Real) const {
	return 0.;
}
Real TmdGaussianSet::E_tilde(Hadron, unsigned, Real, Real) const {
	return 0.;
}

Real TmdGaussianSet::xf1(unsigned fl, Real x, Real Q_sq, Real k_perp_sq) const {
	return xf1(fl, x, Q_sq)*gaussian(k_perp_sq, mean_f1);
}
Real TmdGaussianSet::xf1Tperp(unsigned fl, Real x, Real Q_sq, Real k_perp_sq) const {
	return xf1Tperp(fl, x, Q_sq)*gaussian(k_perp_sq, mean_f1Tperp);
}
Real TmdGaussianSet::xfT(unsigned fl, Real x, Real Q_sq, Real k_perp_sq) const {
	return xfT(fl, x, Q_sq)*gaussian(k_perp_sq, mean_fT);
}
Real TmdGaussianSet::xfperp(unsigned fl, Real x, Real Q_sq, Real k_perp_sq) const {
	return xfperp(fl, x, Q_sq)*gaussian(k_perp_sq, mean_fperp);
}
Real TmdGaussianSet::xfLperp(unsigned fl, Real x, Real Q_sq, Real k_perp_sq) const {
	return xfLperp(fl, x, Q_sq)*gaussian(k_perp_sq, mean_fLperp);
}
Real TmdGaussianSet::xfTperp(unsigned fl, Real x, Real Q_sq, Real k_perp_sq) const {
	return xfTperp(fl, x, Q_sq)*gaussian(k_perp_sq, mean_fTperp);
}
Real TmdGaussianSet::xg1(unsigned fl, Real x, Real Q_sq, Real k_perp_sq) const {
	return xg1(fl, x, Q_sq)*gaussian(k_perp_sq, mean_g1);
}
Real TmdGaussianSet::xg1Tperp(unsigned fl, Real x, Real Q_sq, Real k_perp_sq) const {
	return xg1Tperp(fl, x, Q_sq)*gaussian(k_perp_sq, mean_g1Tperp);
}
Real TmdGaussianSet::xgT(unsigned fl, Real x, Real Q_sq, Real k_perp_sq) const {
	return xgT(fl, x, Q_sq)*gaussian(k_perp_sq, mean_gT);
}
Real TmdGaussianSet::xgperp(unsigned fl, Real x, Real Q_sq, Real k_perp_sq) const {
	return xgperp(fl, x, Q_sq)*gaussian(k_perp_sq, mean_gperp);
}
Real TmdGaussianSet::xgLperp(unsigned fl, Real x, Real Q_sq, Real k_perp_sq) const {
	return xgLperp(fl, x, Q_sq)*gaussian(k_perp_sq, mean_gLperp);
}
Real TmdGaussianSet::xgTperp(unsigned fl, Real x, Real Q_sq, Real k_perp_sq) const {
	return xgTperp(fl, x, Q_sq)*gaussian(k_perp_sq, mean_gTperp);
}
Real TmdGaussianSet::xh1(unsigned fl, Real x, Real Q_sq, Real k_perp_sq) const {
	return xh1(fl, x, Q_sq)*gaussian(k_perp_sq, mean_h1);
}
Real TmdGaussianSet::xh1perp(unsigned fl, Real x, Real Q_sq, Real k_perp_sq) const {
	return xh1perp(fl, x, Q_sq)*gaussian(k_perp_sq, mean_h1perp);
}
Real TmdGaussianSet::xh1Lperp(unsigned fl, Real x, Real Q_sq, Real k_perp_sq) const {
	return xh1Lperp(fl, x, Q_sq)*gaussian(k_perp_sq, mean_h1Lperp);
}
Real TmdGaussianSet::xh1Tperp(unsigned fl, Real x, Real Q_sq, Real k_perp_sq) const {
	return xh1Tperp(fl, x, Q_sq)*gaussian(k_perp_sq, mean_h1Tperp);
}
Real TmdGaussianSet::xh(unsigned fl, Real x, Real Q_sq, Real k_perp_sq) const {
	return xh(fl, x, Q_sq)*gaussian(k_perp_sq, mean_h);
}
Real TmdGaussianSet::xhL(unsigned fl, Real x, Real Q_sq, Real k_perp_sq) const {
	return xhL(fl, x, Q_sq)*gaussian(k_perp_sq, mean_hL);
}
Real TmdGaussianSet::xhT(unsigned fl, Real x, Real Q_sq, Real k_perp_sq) const {
	return xhT(fl, x, Q_sq)*gaussian(k_perp_sq, mean_hT);
}
Real TmdGaussianSet::xhTperp(unsigned fl, Real x, Real Q_sq, Real k_perp_sq) const {
	return xhTperp(fl, x, Q_sq)*gaussian(k_perp_sq, mean_hTperp);
}
Real TmdGaussianSet::xe(unsigned fl, Real x, Real Q_sq, Real k_perp_sq) const {
	return xe(fl, x, Q_sq)*gaussian(k_perp_sq, mean_e);
}
Real TmdGaussianSet::xeL(unsigned fl, Real x, Real Q_sq, Real k_perp_sq) const {
	return xeL(fl, x, Q_sq)*gaussian(k_perp_sq, mean_eL);
}
Real TmdGaussianSet::xeT(unsigned fl, Real x, Real Q_sq, Real k_perp_sq) const {
	return xeT(fl, x, Q_sq)*gaussian(k_perp_sq, mean_eT);
}
Real TmdGaussianSet::xeTperp(unsigned fl, Real x, Real Q_sq, Real k_perp_sq) const {
	return xeTperp(fl, x, Q_sq)*gaussian(k_perp_sq, mean_eTperp);
}

Real TmdGaussianSet::D1(Hadron h, unsigned fl, Real z, Real Q_sq, Real p_perp_sq) const {
	return D1(h, fl, z, Q_sq)*gaussian(p_perp_sq, mean_D1);
}
Real TmdGaussianSet::H1perp(Hadron h, unsigned fl, Real z, Real Q_sq, Real p_perp_sq) const {
	return H1perp(h, fl, z, Q_sq)*gaussian(p_perp_sq, mean_H1perp);
}
Real TmdGaussianSet::Dperp_tilde(Hadron h, unsigned fl, Real z, Real Q_sq, Real p_perp_sq) const {
	return Dperp_tilde(h, fl, z, Q_sq)*gaussian(p_perp_sq, mean_Dperp_tilde);
}
Real TmdGaussianSet::H_tilde(Hadron h, unsigned fl, Real z, Real Q_sq, Real p_perp_sq) const {
	return H_tilde(h, fl, z, Q_sq)*gaussian(p_perp_sq, mean_H_tilde);
}
Real TmdGaussianSet::Gperp_tilde(Hadron h, unsigned fl, Real z, Real Q_sq, Real p_perp_sq) const {
	return Gperp_tilde(h, fl, z, Q_sq)*gaussian(p_perp_sq, mean_Gperp_tilde);
}
Real TmdGaussianSet::E_tilde(Hadron h, unsigned fl, Real z, Real Q_sq, Real p_perp_sq) const {
	return E_tilde(h, fl, z, Q_sq)*gaussian(p_perp_sq, mean_E_tilde);
}

// WW-type approximation.
Real TmdWwSet::xf1(unsigned, Real, Real, Real) const {
	return 0.;
}
Real TmdWwSet::xf1Tperp(unsigned, Real, Real, Real) const {
	return 0.;
}
Real TmdWwSet::xg1(unsigned, Real, Real, Real) const {
	return 0.;
}
Real TmdWwSet::xg1Tperp(unsigned, Real, Real, Real) const {
	return 0.;
}
Real TmdWwSet::xh1(unsigned, Real, Real, Real) const {
	return 0.;
}
Real TmdWwSet::xh1perp(unsigned, Real, Real, Real) const {
	return 0.;
}
Real TmdWwSet::xh1Lperp(unsigned, Real, Real, Real) const {
	return 0.;
}
Real TmdWwSet::xh1Tperp(unsigned, Real, Real, Real) const {
	return 0.;
}

Real TmdWwSet::D1(Hadron, unsigned, Real, Real, Real) const {
	return 0.;
}
Real TmdWwSet::H1perp(Hadron, unsigned, Real, Real, Real) const {
	return 0.;
}

Real TmdWwSet::xf1TperpM1(unsigned fl, Real x, Real Q_sq, Real k_perp_sq) const {
	Real M = mass(target);
	return k_perp_sq/(2.*M*M)*xf1Tperp(fl, x, Q_sq, k_perp_sq);
}
Real TmdWwSet::xg1TperpM1(unsigned fl, Real x, Real Q_sq, Real k_perp_sq) const {
	Real M = mass(target);
	return k_perp_sq/(2.*M*M)*xg1Tperp(fl, x, Q_sq, k_perp_sq);
}
Real TmdWwSet::xh1perpM1(unsigned fl, Real x, Real Q_sq, Real k_perp_sq) const {
	Real M = mass(target);
	return k_perp_sq/(2.*M*M)*xh1perp(fl, x, Q_sq, k_perp_sq);
}
Real TmdWwSet::xh1LperpM1(unsigned fl, Real x, Real Q_sq, Real k_perp_sq) const {
	Real M = mass(target);
	return k_perp_sq/(2.*M*M)*xh1Lperp(fl, x, Q_sq, k_perp_sq);
}
Real TmdWwSet::xh1TperpM1(unsigned fl, Real x, Real Q_sq, Real k_perp_sq) const {
	Real M = mass(target);
	return k_perp_sq/(2.*M*M)*xh1Tperp(fl, x, Q_sq, k_perp_sq);
}

Real TmdWwSet::xfT(unsigned fl, Real x, Real Q_sq, Real k_perp_sq) const {
	return -xf1TperpM1(fl, x, Q_sq, k_perp_sq)/x;
}
Real TmdWwSet::xfperp(unsigned fl, Real x, Real Q_sq, Real k_perp_sq) const {
	return xf1(fl, x, Q_sq, k_perp_sq)/x;
}
Real TmdWwSet::xfLperp(unsigned, Real, Real, Real) const {
	return 0.;
}
Real TmdWwSet::xfTperp(unsigned fl, Real x, Real Q_sq, Real k_perp_sq) const {
	return xf1Tperp(fl, x, Q_sq, k_perp_sq)/x;
}
Real TmdWwSet::xgT(unsigned fl, Real x, Real Q_sq, Real k_perp_sq) const {
	return xg1TperpM1(fl, x, Q_sq, k_perp_sq)/x;
}
Real TmdWwSet::xgperp(unsigned, Real, Real, Real) const {
	return 0.;
}
Real TmdWwSet::xgLperp(unsigned fl, Real x, Real Q_sq, Real k_perp_sq) const {
	return xg1(fl, x, Q_sq, k_perp_sq)/x;
}
Real TmdWwSet::xgTperp(unsigned fl, Real x, Real Q_sq, Real k_perp_sq) const {
	return xg1Tperp(fl, x, Q_sq, k_perp_sq)/x;
}
Real TmdWwSet::xh(unsigned fl, Real x, Real Q_sq, Real k_perp_sq) const {
	return -2.*xh1perpM1(fl, x, Q_sq, k_perp_sq)/x;
}
Real TmdWwSet::xhL(unsigned fl, Real x, Real Q_sq, Real k_perp_sq) const {
	return -2.*xh1LperpM1(fl, x, Q_sq, k_perp_sq)/x;
}
Real TmdWwSet::xhT(unsigned fl, Real x, Real Q_sq, Real k_perp_sq) const {
	return -(xh1(fl, x, Q_sq, k_perp_sq) + xh1TperpM1(fl, x, Q_sq, k_perp_sq))/x;
}
Real TmdWwSet::xhTperp(unsigned fl, Real x, Real Q_sq, Real k_perp_sq) const {
	return (xh1(fl, x, Q_sq, k_perp_sq) - xh1TperpM1(fl, x, Q_sq, k_perp_sq))/x;
}
Real TmdWwSet::xe(unsigned, Real, Real, Real) const {
	return 0.;
}
Real TmdWwSet::xeL(unsigned, Real, Real, Real) const {
	return 0.;
}
Real TmdWwSet::xeT(unsigned, Real, Real, Real) const {
	return 0.;
}
Real TmdWwSet::xeTperp(unsigned, Real, Real, Real) const {
	return 0.;
}

Real TmdWwSet::Dperp_tilde(Hadron, unsigned, Real, Real, Real) const {
	return 0.;
}
Real TmdWwSet::H_tilde(Hadron, unsigned, Real, Real, Real) const {
	return 0.;
}
Real TmdWwSet::Gperp_tilde(Hadron, unsigned, Real, Real, Real) const {
	return 0.;
}
Real TmdWwSet::E_tilde(Hadron, unsigned, Real, Real, Real) const {
	return 0.;
}

// Gaussian and WW-type approximation combined.
TmdGaussianWwSet::TmdGaussianWwSet(
	unsigned flavor_count,
	constant::Nucleus target,
	Real mean_f1,
	Real mean_f1Tperp,
	Real mean_fT,
	Real mean_g1,
	Real mean_g1Tperp,
	Real mean_gT,
	Real mean_h1,
	Real mean_h1perp,
	Real mean_h1Lperp,
	Real mean_h1Tperp,
	Real mean_h,
	Real mean_hL,
	Real mean_hT,
	Real mean_hTperp,
	Real mean_D1,
	Real mean_H1perp) :
	// We choose infinity for the width of the structure functions that are to
	// be neglected.
	TmdGaussianSet(
		flavor_count,
		target,
		mean_f1,
		mean_f1Tperp,
		mean_fT,
		// `mean_fperp`
		mean_f1,
		// `mean_fLperp`
		INF,
		// `mean_fTperp`
		mean_f1Tperp,
		mean_g1,
		mean_g1Tperp,
		mean_gT,
		// `mean_gperp`
		INF,
		// `mean_gLperp`
		mean_g1,
		// `mean_gTperp`
		mean_g1Tperp,
		mean_h1,
		mean_h1perp,
		mean_h1Lperp,
		mean_h1Tperp,
		mean_h,
		mean_hL,
		mean_hT,
		mean_hTperp,
		// `mean_e`
		INF,
		// `mean_eL`
		INF,
		// `mean_eT`
		INF,
		// `mean_eTperp`
		INF,
		mean_D1,
		mean_H1perp,
		// `mean_Dperp_tilde`
		INF,
		// `mean_H_tilde`
		INF,
		// `mean_Gperp_tilde`
		INF,
		// `mean_E_tilde`
		INF) { }

Real TmdGaussianWwSet::xf1(unsigned, Real, Real) const {
	return 0.;
}
Real TmdGaussianWwSet::xf1Tperp(unsigned, Real, Real) const {
	return 0.;
}
Real TmdGaussianWwSet::xg1(unsigned, Real, Real) const {
	return 0.;
}
Real TmdGaussianWwSet::xg1Tperp(unsigned, Real, Real) const {
	return 0.;
}
Real TmdGaussianWwSet::xh1(unsigned, Real, Real) const {
	return 0.;
}
Real TmdGaussianWwSet::xh1perp(unsigned, Real, Real) const {
	return 0.;
}
Real TmdGaussianWwSet::xh1Lperp(unsigned, Real, Real) const {
	return 0.;
}
Real TmdGaussianWwSet::xh1Tperp(unsigned, Real, Real) const {
	return 0.;
}

Real TmdGaussianWwSet::D1(Hadron, unsigned, Real, Real) const {
	return 0.;
}
Real TmdGaussianWwSet::H1perp(Hadron, unsigned, Real, Real) const {
	return 0.;
}

// Following [2], we take the approach of first integrating over `k_perp_sq`
// with the Gaussian approximation, and then using Ww-type approximation. This
// is done to avoid contradictions between these two approximations (see section
// 4.4 of [2]).
Real TmdGaussianWwSet::xf1TperpM1(unsigned fl, Real x, Real Q_sq) const {
	Real M = mass(target);
	return mean_f1Tperp/(2.*M*M)*xf1Tperp(fl, x, Q_sq);
}
Real TmdGaussianWwSet::xg1TperpM1(unsigned fl, Real x, Real Q_sq) const {
	Real M = mass(target);
	return mean_g1Tperp/(2.*M*M)*xg1Tperp(fl, x, Q_sq);
}
Real TmdGaussianWwSet::xh1perpM1(unsigned fl, Real x, Real Q_sq) const {
	Real M = mass(target);
	return mean_h1perp/(2.*M*M)*xh1perp(fl, x, Q_sq);
}
Real TmdGaussianWwSet::xh1LperpM1(unsigned fl, Real x, Real Q_sq) const {
	Real M = mass(target);
	return mean_h1Lperp/(2.*M*M)*xh1Lperp(fl, x, Q_sq);
}
Real TmdGaussianWwSet::xh1TperpM1(unsigned fl, Real x, Real Q_sq) const {
	Real M = mass(target);
	return mean_h1Tperp/(2.*M*M)*xh1Tperp(fl, x, Q_sq);
}

Real TmdGaussianWwSet::xfT(unsigned fl, Real x, Real Q_sq) const {
	return -xf1TperpM1(fl, x, Q_sq)/x;
}
Real TmdGaussianWwSet::xfperp(unsigned fl, Real x, Real Q_sq) const {
	return xf1(fl, x, Q_sq)/x;
}
Real TmdGaussianWwSet::xfLperp(unsigned, Real, Real) const {
	return 0.;
}
Real TmdGaussianWwSet::xfTperp(unsigned fl, Real x, Real Q_sq) const {
	return xf1Tperp(fl, x, Q_sq)/x;
}
Real TmdGaussianWwSet::xgT(unsigned fl, Real x, Real Q_sq) const {
	return xg1TperpM1(fl, x, Q_sq)/x;
}
Real TmdGaussianWwSet::xgperp(unsigned, Real, Real) const {
	return 0.;
}
Real TmdGaussianWwSet::xgLperp(unsigned fl, Real x, Real Q_sq) const {
	return xg1(fl, x, Q_sq)/x;
}
Real TmdGaussianWwSet::xgTperp(unsigned fl, Real x, Real Q_sq) const {
	return xg1Tperp(fl, x, Q_sq)/x;
}
Real TmdGaussianWwSet::xh(unsigned fl, Real x, Real Q_sq) const {
	return -2.*xh1perpM1(fl, x, Q_sq)/x;
}
Real TmdGaussianWwSet::xhL(unsigned fl, Real x, Real Q_sq) const {
	return -2.*xh1LperpM1(fl, x, Q_sq)/x;
}
Real TmdGaussianWwSet::xhT(unsigned fl, Real x, Real Q_sq) const {
	return -(xh1(fl, x, Q_sq) + xh1TperpM1(fl, x, Q_sq))/x;
}
Real TmdGaussianWwSet::xhTperp(unsigned fl, Real x, Real Q_sq) const {
	return (xh1(fl, x, Q_sq) - xh1TperpM1(fl, x, Q_sq))/x;
}
Real TmdGaussianWwSet::xe(unsigned, Real, Real) const {
	return 0.;
}
Real TmdGaussianWwSet::xeL(unsigned, Real, Real) const {
	return 0.;
}
Real TmdGaussianWwSet::xeT(unsigned, Real, Real) const {
	return 0.;
}
Real TmdGaussianWwSet::xeTperp(unsigned, Real, Real) const {
	return 0.;
}

Real TmdGaussianWwSet::Dperp_tilde(Hadron, unsigned, Real, Real) const {
	return 0.;
}
Real TmdGaussianWwSet::H_tilde(Hadron, unsigned, Real, Real) const {
	return 0.;
}
Real TmdGaussianWwSet::Gperp_tilde(Hadron, unsigned, Real, Real) const {
	return 0.;
}
Real TmdGaussianWwSet::E_tilde(Hadron, unsigned, Real, Real) const {
	return 0.;
}

