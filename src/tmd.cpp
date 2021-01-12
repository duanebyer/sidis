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
GaussianTmdSet::GaussianTmdSet(
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

Real GaussianTmdSet::xf1(unsigned, Real, Real) const {
	return 0.;
}
Real GaussianTmdSet::xf1Tperp(unsigned, Real, Real) const {
	return 0.;
}
Real GaussianTmdSet::xfT(unsigned, Real, Real) const {
	return 0.;
}
Real GaussianTmdSet::xfperp(unsigned, Real, Real) const {
	return 0.;
}
Real GaussianTmdSet::xfLperp(unsigned, Real, Real) const {
	return 0.;
}
Real GaussianTmdSet::xfTperp(unsigned, Real, Real) const {
	return 0.;
}
Real GaussianTmdSet::xg1(unsigned, Real, Real) const {
	return 0.;
}
Real GaussianTmdSet::xg1Tperp(unsigned, Real, Real) const {
	return 0.;
}
Real GaussianTmdSet::xgT(unsigned, Real, Real) const {
	return 0.;
}
Real GaussianTmdSet::xgperp(unsigned, Real, Real) const {
	return 0.;
}
Real GaussianTmdSet::xgLperp(unsigned, Real, Real) const {
	return 0.;
}
Real GaussianTmdSet::xgTperp(unsigned, Real, Real) const {
	return 0.;
}
Real GaussianTmdSet::xh1(unsigned, Real, Real) const {
	return 0.;
}
Real GaussianTmdSet::xh1perp(unsigned, Real, Real) const {
	return 0.;
}
Real GaussianTmdSet::xh1Lperp(unsigned, Real, Real) const {
	return 0.;
}
Real GaussianTmdSet::xh1Tperp(unsigned, Real, Real) const {
	return 0.;
}
Real GaussianTmdSet::xh(unsigned, Real, Real) const {
	return 0.;
}
Real GaussianTmdSet::xhL(unsigned, Real, Real) const {
	return 0.;
}
Real GaussianTmdSet::xhT(unsigned, Real, Real) const {
	return 0.;
}
Real GaussianTmdSet::xhTperp(unsigned, Real, Real) const {
	return 0.;
}
Real GaussianTmdSet::xe(unsigned, Real, Real) const {
	return 0.;
}
Real GaussianTmdSet::xeL(unsigned, Real, Real) const {
	return 0.;
}
Real GaussianTmdSet::xeT(unsigned, Real, Real) const {
	return 0.;
}
Real GaussianTmdSet::xeTperp(unsigned, Real, Real) const {
	return 0.;
}

Real GaussianTmdSet::D1(Hadron, unsigned, Real, Real) const {
	return 0.;
}
Real GaussianTmdSet::H1perp(Hadron, unsigned, Real, Real) const {
	return 0.;
}
Real GaussianTmdSet::Dperp_tilde(Hadron, unsigned, Real, Real) const {
	return 0.;
}
Real GaussianTmdSet::H_tilde(Hadron, unsigned, Real, Real) const {
	return 0.;
}
Real GaussianTmdSet::Gperp_tilde(Hadron, unsigned, Real, Real) const {
	return 0.;
}
Real GaussianTmdSet::E_tilde(Hadron, unsigned, Real, Real) const {
	return 0.;
}

Real GaussianTmdSet::xf1(unsigned fl, Real x, Real Q_sq, Real k_perp_sq) const {
	return xf1(fl, x, Q_sq)*gaussian(k_perp_sq, mean_f1);
}
Real GaussianTmdSet::xf1Tperp(unsigned fl, Real x, Real Q_sq, Real k_perp_sq) const {
	return xf1Tperp(fl, x, Q_sq)*gaussian(k_perp_sq, mean_f1Tperp);
}
Real GaussianTmdSet::xfT(unsigned fl, Real x, Real Q_sq, Real k_perp_sq) const {
	return xfT(fl, x, Q_sq)*gaussian(k_perp_sq, mean_fT);
}
Real GaussianTmdSet::xfperp(unsigned fl, Real x, Real Q_sq, Real k_perp_sq) const {
	return xfperp(fl, x, Q_sq)*gaussian(k_perp_sq, mean_fperp);
}
Real GaussianTmdSet::xfLperp(unsigned fl, Real x, Real Q_sq, Real k_perp_sq) const {
	return xfLperp(fl, x, Q_sq)*gaussian(k_perp_sq, mean_fLperp);
}
Real GaussianTmdSet::xfTperp(unsigned fl, Real x, Real Q_sq, Real k_perp_sq) const {
	return xfTperp(fl, x, Q_sq)*gaussian(k_perp_sq, mean_fTperp);
}
Real GaussianTmdSet::xg1(unsigned fl, Real x, Real Q_sq, Real k_perp_sq) const {
	return xg1(fl, x, Q_sq)*gaussian(k_perp_sq, mean_g1);
}
Real GaussianTmdSet::xg1Tperp(unsigned fl, Real x, Real Q_sq, Real k_perp_sq) const {
	return xg1Tperp(fl, x, Q_sq)*gaussian(k_perp_sq, mean_g1Tperp);
}
Real GaussianTmdSet::xgT(unsigned fl, Real x, Real Q_sq, Real k_perp_sq) const {
	return xgT(fl, x, Q_sq)*gaussian(k_perp_sq, mean_gT);
}
Real GaussianTmdSet::xgperp(unsigned fl, Real x, Real Q_sq, Real k_perp_sq) const {
	return xgperp(fl, x, Q_sq)*gaussian(k_perp_sq, mean_gperp);
}
Real GaussianTmdSet::xgLperp(unsigned fl, Real x, Real Q_sq, Real k_perp_sq) const {
	return xgLperp(fl, x, Q_sq)*gaussian(k_perp_sq, mean_gLperp);
}
Real GaussianTmdSet::xgTperp(unsigned fl, Real x, Real Q_sq, Real k_perp_sq) const {
	return xgTperp(fl, x, Q_sq)*gaussian(k_perp_sq, mean_gTperp);
}
Real GaussianTmdSet::xh1(unsigned fl, Real x, Real Q_sq, Real k_perp_sq) const {
	return xh1(fl, x, Q_sq)*gaussian(k_perp_sq, mean_h1);
}
Real GaussianTmdSet::xh1perp(unsigned fl, Real x, Real Q_sq, Real k_perp_sq) const {
	return xh1perp(fl, x, Q_sq)*gaussian(k_perp_sq, mean_h1perp);
}
Real GaussianTmdSet::xh1Lperp(unsigned fl, Real x, Real Q_sq, Real k_perp_sq) const {
	return xh1Lperp(fl, x, Q_sq)*gaussian(k_perp_sq, mean_h1Lperp);
}
Real GaussianTmdSet::xh1Tperp(unsigned fl, Real x, Real Q_sq, Real k_perp_sq) const {
	return xh1Tperp(fl, x, Q_sq)*gaussian(k_perp_sq, mean_h1Tperp);
}
Real GaussianTmdSet::xh(unsigned fl, Real x, Real Q_sq, Real k_perp_sq) const {
	return xh(fl, x, Q_sq)*gaussian(k_perp_sq, mean_h);
}
Real GaussianTmdSet::xhL(unsigned fl, Real x, Real Q_sq, Real k_perp_sq) const {
	return xhL(fl, x, Q_sq)*gaussian(k_perp_sq, mean_hL);
}
Real GaussianTmdSet::xhT(unsigned fl, Real x, Real Q_sq, Real k_perp_sq) const {
	return xhT(fl, x, Q_sq)*gaussian(k_perp_sq, mean_hT);
}
Real GaussianTmdSet::xhTperp(unsigned fl, Real x, Real Q_sq, Real k_perp_sq) const {
	return xhTperp(fl, x, Q_sq)*gaussian(k_perp_sq, mean_hTperp);
}
Real GaussianTmdSet::xe(unsigned fl, Real x, Real Q_sq, Real k_perp_sq) const {
	return xe(fl, x, Q_sq)*gaussian(k_perp_sq, mean_e);
}
Real GaussianTmdSet::xeL(unsigned fl, Real x, Real Q_sq, Real k_perp_sq) const {
	return xeL(fl, x, Q_sq)*gaussian(k_perp_sq, mean_eL);
}
Real GaussianTmdSet::xeT(unsigned fl, Real x, Real Q_sq, Real k_perp_sq) const {
	return xeT(fl, x, Q_sq)*gaussian(k_perp_sq, mean_eT);
}
Real GaussianTmdSet::xeTperp(unsigned fl, Real x, Real Q_sq, Real k_perp_sq) const {
	return xeTperp(fl, x, Q_sq)*gaussian(k_perp_sq, mean_eTperp);
}

Real GaussianTmdSet::D1(Hadron h, unsigned fl, Real z, Real Q_sq, Real p_perp_sq) const {
	return D1(h, fl, z, Q_sq)*gaussian(p_perp_sq, mean_D1);
}
Real GaussianTmdSet::H1perp(Hadron h, unsigned fl, Real z, Real Q_sq, Real p_perp_sq) const {
	return H1perp(h, fl, z, Q_sq)*gaussian(p_perp_sq, mean_H1perp);
}
Real GaussianTmdSet::Dperp_tilde(Hadron h, unsigned fl, Real z, Real Q_sq, Real p_perp_sq) const {
	return Dperp_tilde(h, fl, z, Q_sq)*gaussian(p_perp_sq, mean_Dperp_tilde);
}
Real GaussianTmdSet::H_tilde(Hadron h, unsigned fl, Real z, Real Q_sq, Real p_perp_sq) const {
	return H_tilde(h, fl, z, Q_sq)*gaussian(p_perp_sq, mean_H_tilde);
}
Real GaussianTmdSet::Gperp_tilde(Hadron h, unsigned fl, Real z, Real Q_sq, Real p_perp_sq) const {
	return Gperp_tilde(h, fl, z, Q_sq)*gaussian(p_perp_sq, mean_Gperp_tilde);
}
Real GaussianTmdSet::E_tilde(Hadron h, unsigned fl, Real z, Real Q_sq, Real p_perp_sq) const {
	return E_tilde(h, fl, z, Q_sq)*gaussian(p_perp_sq, mean_E_tilde);
}

// WW-type approximation.
Real WwTmdSet::xf1(unsigned, Real, Real, Real) const {
	return 0.;
}
Real WwTmdSet::xf1Tperp(unsigned, Real, Real, Real) const {
	return 0.;
}
Real WwTmdSet::xg1(unsigned, Real, Real, Real) const {
	return 0.;
}
Real WwTmdSet::xg1Tperp(unsigned, Real, Real, Real) const {
	return 0.;
}
Real WwTmdSet::xh1(unsigned, Real, Real, Real) const {
	return 0.;
}
Real WwTmdSet::xh1perp(unsigned, Real, Real, Real) const {
	return 0.;
}
Real WwTmdSet::xh1Lperp(unsigned, Real, Real, Real) const {
	return 0.;
}
Real WwTmdSet::xh1Tperp(unsigned, Real, Real, Real) const {
	return 0.;
}

Real WwTmdSet::D1(Hadron, unsigned, Real, Real, Real) const {
	return 0.;
}
Real WwTmdSet::H1perp(Hadron, unsigned, Real, Real, Real) const {
	return 0.;
}

Real WwTmdSet::xf1TperpM1(unsigned fl, Real x, Real Q_sq, Real k_perp_sq) const {
	Real M = mass(target);
	return k_perp_sq/(2.*M*M)*xf1Tperp(fl, x, Q_sq, k_perp_sq);
}
Real WwTmdSet::xg1TperpM1(unsigned fl, Real x, Real Q_sq, Real k_perp_sq) const {
	Real M = mass(target);
	return k_perp_sq/(2.*M*M)*xg1Tperp(fl, x, Q_sq, k_perp_sq);
}
Real WwTmdSet::xh1perpM1(unsigned fl, Real x, Real Q_sq, Real k_perp_sq) const {
	Real M = mass(target);
	return k_perp_sq/(2.*M*M)*xh1perp(fl, x, Q_sq, k_perp_sq);
}
Real WwTmdSet::xh1LperpM1(unsigned fl, Real x, Real Q_sq, Real k_perp_sq) const {
	Real M = mass(target);
	return k_perp_sq/(2.*M*M)*xh1Lperp(fl, x, Q_sq, k_perp_sq);
}
Real WwTmdSet::xh1TperpM1(unsigned fl, Real x, Real Q_sq, Real k_perp_sq) const {
	Real M = mass(target);
	return k_perp_sq/(2.*M*M)*xh1Tperp(fl, x, Q_sq, k_perp_sq);
}

Real WwTmdSet::xfT(unsigned fl, Real x, Real Q_sq, Real k_perp_sq) const {
	return -xf1TperpM1(fl, x, Q_sq, k_perp_sq)/x;
}
Real WwTmdSet::xfperp(unsigned fl, Real x, Real Q_sq, Real k_perp_sq) const {
	return xf1(fl, x, Q_sq, k_perp_sq)/x;
}
Real WwTmdSet::xfLperp(unsigned, Real, Real, Real) const {
	return 0.;
}
Real WwTmdSet::xfTperp(unsigned fl, Real x, Real Q_sq, Real k_perp_sq) const {
	return xf1Tperp(fl, x, Q_sq, k_perp_sq)/x;
}
Real WwTmdSet::xgT(unsigned fl, Real x, Real Q_sq, Real k_perp_sq) const {
	return xg1TperpM1(fl, x, Q_sq, k_perp_sq)/x;
}
Real WwTmdSet::xgperp(unsigned, Real, Real, Real) const {
	return 0.;
}
Real WwTmdSet::xgLperp(unsigned fl, Real x, Real Q_sq, Real k_perp_sq) const {
	return xg1(fl, x, Q_sq, k_perp_sq)/x;
}
Real WwTmdSet::xgTperp(unsigned fl, Real x, Real Q_sq, Real k_perp_sq) const {
	return xg1Tperp(fl, x, Q_sq, k_perp_sq)/x;
}
Real WwTmdSet::xh(unsigned fl, Real x, Real Q_sq, Real k_perp_sq) const {
	return -2.*xh1perpM1(fl, x, Q_sq, k_perp_sq)/x;
}
Real WwTmdSet::xhL(unsigned fl, Real x, Real Q_sq, Real k_perp_sq) const {
	return -2.*xh1LperpM1(fl, x, Q_sq, k_perp_sq)/x;
}
Real WwTmdSet::xhT(unsigned fl, Real x, Real Q_sq, Real k_perp_sq) const {
	return -(xh1(fl, x, Q_sq, k_perp_sq) + xh1TperpM1(fl, x, Q_sq, k_perp_sq))/x;
}
Real WwTmdSet::xhTperp(unsigned fl, Real x, Real Q_sq, Real k_perp_sq) const {
	return (xh1(fl, x, Q_sq, k_perp_sq) - xh1TperpM1(fl, x, Q_sq, k_perp_sq))/x;
}
Real WwTmdSet::xe(unsigned, Real, Real, Real) const {
	return 0.;
}
Real WwTmdSet::xeL(unsigned, Real, Real, Real) const {
	return 0.;
}
Real WwTmdSet::xeT(unsigned, Real, Real, Real) const {
	return 0.;
}
Real WwTmdSet::xeTperp(unsigned, Real, Real, Real) const {
	return 0.;
}

Real WwTmdSet::Dperp_tilde(Hadron, unsigned, Real, Real, Real) const {
	return 0.;
}
Real WwTmdSet::H_tilde(Hadron, unsigned, Real, Real, Real) const {
	return 0.;
}
Real WwTmdSet::Gperp_tilde(Hadron, unsigned, Real, Real, Real) const {
	return 0.;
}
Real WwTmdSet::E_tilde(Hadron, unsigned, Real, Real, Real) const {
	return 0.;
}

// Gaussian and WW-type approximation combined.
GaussianWwTmdSet::GaussianWwTmdSet(
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
	GaussianTmdSet(
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

Real GaussianWwTmdSet::xf1(unsigned, Real, Real) const {
	return 0.;
}
Real GaussianWwTmdSet::xf1Tperp(unsigned, Real, Real) const {
	return 0.;
}
Real GaussianWwTmdSet::xg1(unsigned, Real, Real) const {
	return 0.;
}
Real GaussianWwTmdSet::xg1Tperp(unsigned, Real, Real) const {
	return 0.;
}
Real GaussianWwTmdSet::xh1(unsigned, Real, Real) const {
	return 0.;
}
Real GaussianWwTmdSet::xh1perp(unsigned, Real, Real) const {
	return 0.;
}
Real GaussianWwTmdSet::xh1Lperp(unsigned, Real, Real) const {
	return 0.;
}
Real GaussianWwTmdSet::xh1Tperp(unsigned, Real, Real) const {
	return 0.;
}

Real GaussianWwTmdSet::D1(Hadron, unsigned, Real, Real) const {
	return 0.;
}
Real GaussianWwTmdSet::H1perp(Hadron, unsigned, Real, Real) const {
	return 0.;
}

// Following [2], we take the approach of first integrating over `k_perp_sq`
// with the Gaussian approximation, and then using Ww-type approximation. This
// is done to avoid contradictions between these two approximations (see section
// 4.4 of [2]).
Real GaussianWwTmdSet::xf1TperpM1(unsigned fl, Real x, Real Q_sq) const {
	Real M = mass(target);
	return mean_f1Tperp/(2.*M*M)*xf1Tperp(fl, x, Q_sq);
}
Real GaussianWwTmdSet::xg1TperpM1(unsigned fl, Real x, Real Q_sq) const {
	Real M = mass(target);
	return mean_g1Tperp/(2.*M*M)*xg1Tperp(fl, x, Q_sq);
}
Real GaussianWwTmdSet::xh1perpM1(unsigned fl, Real x, Real Q_sq) const {
	Real M = mass(target);
	return mean_h1perp/(2.*M*M)*xh1perp(fl, x, Q_sq);
}
Real GaussianWwTmdSet::xh1LperpM1(unsigned fl, Real x, Real Q_sq) const {
	Real M = mass(target);
	return mean_h1Lperp/(2.*M*M)*xh1Lperp(fl, x, Q_sq);
}
Real GaussianWwTmdSet::xh1TperpM1(unsigned fl, Real x, Real Q_sq) const {
	Real M = mass(target);
	return mean_h1Tperp/(2.*M*M)*xh1Tperp(fl, x, Q_sq);
}

Real GaussianWwTmdSet::xfT(unsigned fl, Real x, Real Q_sq) const {
	return -xf1TperpM1(fl, x, Q_sq)/x;
}
Real GaussianWwTmdSet::xfperp(unsigned fl, Real x, Real Q_sq) const {
	return xf1(fl, x, Q_sq)/x;
}
Real GaussianWwTmdSet::xfLperp(unsigned, Real, Real) const {
	return 0.;
}
Real GaussianWwTmdSet::xfTperp(unsigned fl, Real x, Real Q_sq) const {
	return xf1Tperp(fl, x, Q_sq)/x;
}
Real GaussianWwTmdSet::xgT(unsigned fl, Real x, Real Q_sq) const {
	return xg1TperpM1(fl, x, Q_sq)/x;
}
Real GaussianWwTmdSet::xgperp(unsigned, Real, Real) const {
	return 0.;
}
Real GaussianWwTmdSet::xgLperp(unsigned fl, Real x, Real Q_sq) const {
	return xg1(fl, x, Q_sq)/x;
}
Real GaussianWwTmdSet::xgTperp(unsigned fl, Real x, Real Q_sq) const {
	return xg1Tperp(fl, x, Q_sq)/x;
}
Real GaussianWwTmdSet::xh(unsigned fl, Real x, Real Q_sq) const {
	return -2.*xh1perpM1(fl, x, Q_sq)/x;
}
Real GaussianWwTmdSet::xhL(unsigned fl, Real x, Real Q_sq) const {
	return -2.*xh1LperpM1(fl, x, Q_sq)/x;
}
Real GaussianWwTmdSet::xhT(unsigned fl, Real x, Real Q_sq) const {
	return -(xh1(fl, x, Q_sq) + xh1TperpM1(fl, x, Q_sq))/x;
}
Real GaussianWwTmdSet::xhTperp(unsigned fl, Real x, Real Q_sq) const {
	return (xh1(fl, x, Q_sq) - xh1TperpM1(fl, x, Q_sq))/x;
}
Real GaussianWwTmdSet::xe(unsigned, Real, Real) const {
	return 0.;
}
Real GaussianWwTmdSet::xeL(unsigned, Real, Real) const {
	return 0.;
}
Real GaussianWwTmdSet::xeT(unsigned, Real, Real) const {
	return 0.;
}
Real GaussianWwTmdSet::xeTperp(unsigned, Real, Real) const {
	return 0.;
}

Real GaussianWwTmdSet::Dperp_tilde(Hadron, unsigned, Real, Real) const {
	return 0.;
}
Real GaussianWwTmdSet::H_tilde(Hadron, unsigned, Real, Real) const {
	return 0.;
}
Real GaussianWwTmdSet::Gperp_tilde(Hadron, unsigned, Real, Real) const {
	return 0.;
}
Real GaussianWwTmdSet::E_tilde(Hadron, unsigned, Real, Real) const {
	return 0.;
}

