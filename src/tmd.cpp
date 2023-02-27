#include "sidis/tmd.hpp"

#include <cmath>

#include "sidis/extra/exception.hpp"

using namespace sidis;
using namespace sidis::sf;

TmdSet::TmdSet(unsigned flavor_count, FlavorVec const& charges, part::Nucleus target) :
		flavor_count(flavor_count),
		charges(charges),
		target(target) {
	if (!(flavor_count < MAX_FLAVOR_VEC_SIZE)) {
		throw FlavorOutOfRange(flavor_count);
	}
	if (charges.size() != flavor_count) {
		throw FlavorVecUnexpectedSize(charges.size(), flavor_count);
	}
}

// Standard TMDs.
FlavorVec TmdSet::xf1(Real, Real, Real) const {
	return FlavorVec(flavor_count);
}
FlavorVec TmdSet::xf1Tperp(Real, Real, Real) const {
	return FlavorVec(flavor_count);
}
FlavorVec TmdSet::xfT(Real, Real, Real) const {
	return FlavorVec(flavor_count);
}
FlavorVec TmdSet::xfperp(Real, Real, Real) const {
	return FlavorVec(flavor_count);
}
FlavorVec TmdSet::xfLperp(Real, Real, Real) const {
	return FlavorVec(flavor_count);
}
FlavorVec TmdSet::xfTperp(Real, Real, Real) const {
	return FlavorVec(flavor_count);
}
FlavorVec TmdSet::xg1(Real, Real, Real) const {
	return FlavorVec(flavor_count);
}
FlavorVec TmdSet::xg1Tperp(Real, Real, Real) const {
	return FlavorVec(flavor_count);
}
FlavorVec TmdSet::xgT(Real, Real, Real) const {
	return FlavorVec(flavor_count);
}
FlavorVec TmdSet::xgperp(Real, Real, Real) const {
	return FlavorVec(flavor_count);
}
FlavorVec TmdSet::xgLperp(Real, Real, Real) const {
	return FlavorVec(flavor_count);
}
FlavorVec TmdSet::xgTperp(Real, Real, Real) const {
	return FlavorVec(flavor_count);
}
FlavorVec TmdSet::xh1(Real, Real, Real) const {
	return FlavorVec(flavor_count);
}
FlavorVec TmdSet::xh1perp(Real, Real, Real) const {
	return FlavorVec(flavor_count);
}
FlavorVec TmdSet::xh1Lperp(Real, Real, Real) const {
	return FlavorVec(flavor_count);
}
FlavorVec TmdSet::xh1Tperp(Real, Real, Real) const {
	return FlavorVec(flavor_count);
}
FlavorVec TmdSet::xh(Real, Real, Real) const {
	return FlavorVec(flavor_count);
}
FlavorVec TmdSet::xhL(Real, Real, Real) const {
	return FlavorVec(flavor_count);
}
FlavorVec TmdSet::xhT(Real, Real, Real) const {
	return FlavorVec(flavor_count);
}
FlavorVec TmdSet::xhTperp(Real, Real, Real) const {
	return FlavorVec(flavor_count);
}
FlavorVec TmdSet::xe(Real, Real, Real) const {
	return FlavorVec(flavor_count);
}
FlavorVec TmdSet::xeL(Real, Real, Real) const {
	return FlavorVec(flavor_count);
}
FlavorVec TmdSet::xeT(Real, Real, Real) const {
	return FlavorVec(flavor_count);
}
FlavorVec TmdSet::xeTperp(Real, Real, Real) const {
	return FlavorVec(flavor_count);
}

FlavorVec TmdSet::D1(part::Hadron, Real, Real, Real) const {
	return FlavorVec(flavor_count);
}
FlavorVec TmdSet::H1perp(part::Hadron, Real, Real, Real) const {
	return FlavorVec(flavor_count);
}
FlavorVec TmdSet::Dperp_tilde(part::Hadron, Real, Real, Real) const {
	return FlavorVec(flavor_count);
}
FlavorVec TmdSet::H_tilde(part::Hadron, Real, Real, Real) const {
	return FlavorVec(flavor_count);
}
FlavorVec TmdSet::Gperp_tilde(part::Hadron, Real, Real, Real) const {
	return FlavorVec(flavor_count);
}
FlavorVec TmdSet::E_tilde(part::Hadron, Real, Real, Real) const {
	return FlavorVec(flavor_count);
}

GaussianTmdVars::GaussianTmdVars() :
	f1(1, INF),
	f1Tperp(1, INF),
	fT(1, INF),
	fperp(1, INF),
	fLperp(1, INF),
	fTperp(1, INF),
	g1(1, INF),
	g1Tperp(1, INF),
	gT(1, INF),
	gperp(1, INF),
	gLperp(1, INF),
	gTperp(1, INF),
	h1(1, INF),
	h1perp(1, INF),
	h1Lperp(1, INF),
	h1Tperp(1, INF),
	h(1, INF),
	hL(1, INF),
	hT(1, INF),
	hTperp(1, INF),
	e(1, INF),
	eL(1, INF),
	eT(1, INF),
	eTperp(1, INF),
	D1(1, INF),
	H1perp(1, INF),
	Dperp_tilde(1, INF),
	H_tilde(1, INF),
	Gperp_tilde(1, INF),
	E_tilde(1, INF) { }

GaussianTmdVars::GaussianTmdVars(GaussianWwTmdVars const& ww_vars) :
	// We choose infinity for the width of the structure functions that are to
	// be neglected.
	f1(ww_vars.f1),
	f1Tperp(ww_vars.f1Tperp),
	fT(ww_vars.fT),
	fperp(ww_vars.f1),
	fLperp(f1.size(), INF),
	fTperp(ww_vars.f1Tperp),
	g1(ww_vars.g1),
	g1Tperp(ww_vars.g1Tperp),
	gT(ww_vars.gT),
	gperp(f1.size(), INF),
	gLperp(ww_vars.g1),
	gTperp(ww_vars.g1Tperp),
	h1(ww_vars.h1),
	h1perp(ww_vars.h1perp),
	h1Lperp(ww_vars.h1Lperp),
	h1Tperp(ww_vars.h1Tperp),
	h(ww_vars.h),
	hL(ww_vars.hL),
	hT(ww_vars.hT),
	hTperp(ww_vars.hTperp),
	e(f1.size(), INF),
	eL(f1.size(), INF),
	eT(f1.size(), INF),
	eTperp(f1.size(), INF),
	D1(ww_vars.D1),
	H1perp(ww_vars.H1perp),
	Dperp_tilde(f1.size(), INF),
	H_tilde(f1.size(), INF),
	Gperp_tilde(f1.size(), INF),
	E_tilde(f1.size(), INF) { }

#define FILL_AND_CHECK_TMD_VAR(f) do { \
	if (vars.f.size() == 1) { \
		vars.f = FlavorVec(flavor_count, vars.f[0]); \
	} \
	if (vars.f.size() != flavor_count) { \
		throw FlavorVecUnexpectedSize(vars.f.size(), flavor_count); \
	} \
} while (false)

namespace {
// Checks that all of the variances have one entry for every flavor. If any
// variances have only a single entry, then that entry is copied `flavor_count`
// times. Otherwise, an exception is thrown.
GaussianTmdVars fill_vars(unsigned flavor_count, GaussianTmdVars vars) {
	FILL_AND_CHECK_TMD_VAR(f1);
	FILL_AND_CHECK_TMD_VAR(f1Tperp);
	FILL_AND_CHECK_TMD_VAR(fT);
	FILL_AND_CHECK_TMD_VAR(fperp);
	FILL_AND_CHECK_TMD_VAR(fLperp);
	FILL_AND_CHECK_TMD_VAR(fTperp);
	FILL_AND_CHECK_TMD_VAR(g1);
	FILL_AND_CHECK_TMD_VAR(g1Tperp);
	FILL_AND_CHECK_TMD_VAR(gT);
	FILL_AND_CHECK_TMD_VAR(gperp);
	FILL_AND_CHECK_TMD_VAR(gLperp);
	FILL_AND_CHECK_TMD_VAR(gTperp);
	FILL_AND_CHECK_TMD_VAR(h1);
	FILL_AND_CHECK_TMD_VAR(h1perp);
	FILL_AND_CHECK_TMD_VAR(h1Lperp);
	FILL_AND_CHECK_TMD_VAR(h1Tperp);
	FILL_AND_CHECK_TMD_VAR(h);
	FILL_AND_CHECK_TMD_VAR(hL);
	FILL_AND_CHECK_TMD_VAR(hT);
	FILL_AND_CHECK_TMD_VAR(hTperp);
	FILL_AND_CHECK_TMD_VAR(e);
	FILL_AND_CHECK_TMD_VAR(eL);
	FILL_AND_CHECK_TMD_VAR(eT);
	FILL_AND_CHECK_TMD_VAR(eTperp);
	FILL_AND_CHECK_TMD_VAR(D1);
	FILL_AND_CHECK_TMD_VAR(H1perp);
	FILL_AND_CHECK_TMD_VAR(Dperp_tilde);
	FILL_AND_CHECK_TMD_VAR(H_tilde);
	FILL_AND_CHECK_TMD_VAR(Gperp_tilde);
	FILL_AND_CHECK_TMD_VAR(E_tilde);
	return vars;
}
}

// Gaussian approximation.
GaussianTmdSet::GaussianTmdSet(
		unsigned flavor_count,
		FlavorVec const& charges,
		part::Nucleus target,
		GaussianTmdVars const& vars) :
		TmdSet(flavor_count, charges, target),
		vars(fill_vars(flavor_count, vars)) { }

FlavorVec GaussianTmdSet::xf1(Real, Real) const {
	return FlavorVec(flavor_count);
}
FlavorVec GaussianTmdSet::xf1Tperp(Real, Real) const {
	return FlavorVec(flavor_count);
}
FlavorVec GaussianTmdSet::xfT(Real, Real) const {
	return FlavorVec(flavor_count);
}
FlavorVec GaussianTmdSet::xfperp(Real, Real) const {
	return FlavorVec(flavor_count);
}
FlavorVec GaussianTmdSet::xfLperp(Real, Real) const {
	return FlavorVec(flavor_count);
}
FlavorVec GaussianTmdSet::xfTperp(Real, Real) const {
	return FlavorVec(flavor_count);
}
FlavorVec GaussianTmdSet::xg1(Real, Real) const {
	return FlavorVec(flavor_count);
}
FlavorVec GaussianTmdSet::xg1Tperp(Real, Real) const {
	return FlavorVec(flavor_count);
}
FlavorVec GaussianTmdSet::xgT(Real, Real) const {
	return FlavorVec(flavor_count);
}
FlavorVec GaussianTmdSet::xgperp(Real, Real) const {
	return FlavorVec(flavor_count);
}
FlavorVec GaussianTmdSet::xgLperp(Real, Real) const {
	return FlavorVec(flavor_count);
}
FlavorVec GaussianTmdSet::xgTperp(Real, Real) const {
	return FlavorVec(flavor_count);
}
FlavorVec GaussianTmdSet::xh1(Real, Real) const {
	return FlavorVec(flavor_count);
}
FlavorVec GaussianTmdSet::xh1perp(Real, Real) const {
	return FlavorVec(flavor_count);
}
FlavorVec GaussianTmdSet::xh1Lperp(Real, Real) const {
	return FlavorVec(flavor_count);
}
FlavorVec GaussianTmdSet::xh1Tperp(Real, Real) const {
	return FlavorVec(flavor_count);
}
FlavorVec GaussianTmdSet::xh(Real, Real) const {
	return FlavorVec(flavor_count);
}
FlavorVec GaussianTmdSet::xhL(Real, Real) const {
	return FlavorVec(flavor_count);
}
FlavorVec GaussianTmdSet::xhT(Real, Real) const {
	return FlavorVec(flavor_count);
}
FlavorVec GaussianTmdSet::xhTperp(Real, Real) const {
	return FlavorVec(flavor_count);
}
FlavorVec GaussianTmdSet::xe(Real, Real) const {
	return FlavorVec(flavor_count);
}
FlavorVec GaussianTmdSet::xeL(Real, Real) const {
	return FlavorVec(flavor_count);
}
FlavorVec GaussianTmdSet::xeT(Real, Real) const {
	return FlavorVec(flavor_count);
}
FlavorVec GaussianTmdSet::xeTperp(Real, Real) const {
	return FlavorVec(flavor_count);
}

FlavorVec GaussianTmdSet::D1(part::Hadron, Real, Real) const {
	return FlavorVec(flavor_count);
}
FlavorVec GaussianTmdSet::H1perp(part::Hadron, Real, Real) const {
	return FlavorVec(flavor_count);
}
FlavorVec GaussianTmdSet::Dperp_tilde(part::Hadron, Real, Real) const {
	return FlavorVec(flavor_count);
}
FlavorVec GaussianTmdSet::H_tilde(part::Hadron, Real, Real) const {
	return FlavorVec(flavor_count);
}
FlavorVec GaussianTmdSet::Gperp_tilde(part::Hadron, Real, Real) const {
	return FlavorVec(flavor_count);
}
FlavorVec GaussianTmdSet::E_tilde(part::Hadron, Real, Real) const {
	return FlavorVec(flavor_count);
}

FlavorVec GaussianTmdSet::xf1(Real x, Real Q_sq, Real k_perp_sq) const {
	return tmd_gaussian_factor(vars.f1, k_perp_sq)*xf1(x, Q_sq);
}
FlavorVec GaussianTmdSet::xf1Tperp(Real x, Real Q_sq, Real k_perp_sq) const {
	return tmd_gaussian_factor(vars.f1Tperp, k_perp_sq)*xf1Tperp(x, Q_sq);
}
FlavorVec GaussianTmdSet::xfT(Real x, Real Q_sq, Real k_perp_sq) const {
	return tmd_gaussian_factor(vars.fT, k_perp_sq)*xfT(x, Q_sq);
}
FlavorVec GaussianTmdSet::xfperp(Real x, Real Q_sq, Real k_perp_sq) const {
	return tmd_gaussian_factor(vars.fperp, k_perp_sq)*xfperp(x, Q_sq);
}
FlavorVec GaussianTmdSet::xfLperp(Real x, Real Q_sq, Real k_perp_sq) const {
	return tmd_gaussian_factor(vars.fLperp, k_perp_sq)*xfLperp(x, Q_sq);
}
FlavorVec GaussianTmdSet::xfTperp(Real x, Real Q_sq, Real k_perp_sq) const {
	return tmd_gaussian_factor(vars.fTperp, k_perp_sq)*xfTperp(x, Q_sq);
}
FlavorVec GaussianTmdSet::xg1(Real x, Real Q_sq, Real k_perp_sq) const {
	return tmd_gaussian_factor(vars.g1, k_perp_sq)*xg1(x, Q_sq);
}
FlavorVec GaussianTmdSet::xg1Tperp(Real x, Real Q_sq, Real k_perp_sq) const {
	return tmd_gaussian_factor(vars.g1Tperp, k_perp_sq)*xg1Tperp(x, Q_sq);
}
FlavorVec GaussianTmdSet::xgT(Real x, Real Q_sq, Real k_perp_sq) const {
	return tmd_gaussian_factor(vars.gT, k_perp_sq)*xgT(x, Q_sq);
}
FlavorVec GaussianTmdSet::xgperp(Real x, Real Q_sq, Real k_perp_sq) const {
	return tmd_gaussian_factor(vars.gperp, k_perp_sq)*xgperp(x, Q_sq);
}
FlavorVec GaussianTmdSet::xgLperp(Real x, Real Q_sq, Real k_perp_sq) const {
	return tmd_gaussian_factor(vars.gLperp, k_perp_sq)*xgLperp(x, Q_sq);
}
FlavorVec GaussianTmdSet::xgTperp(Real x, Real Q_sq, Real k_perp_sq) const {
	return tmd_gaussian_factor(vars.gTperp, k_perp_sq)*xgTperp(x, Q_sq);
}
FlavorVec GaussianTmdSet::xh1(Real x, Real Q_sq, Real k_perp_sq) const {
	return tmd_gaussian_factor(vars.h1, k_perp_sq)*xh1(x, Q_sq);
}
FlavorVec GaussianTmdSet::xh1perp(Real x, Real Q_sq, Real k_perp_sq) const {
	return tmd_gaussian_factor(vars.h1perp, k_perp_sq)*xh1perp(x, Q_sq);
}
FlavorVec GaussianTmdSet::xh1Lperp(Real x, Real Q_sq, Real k_perp_sq) const {
	return tmd_gaussian_factor(vars.h1Lperp, k_perp_sq)*xh1Lperp(x, Q_sq);
}
FlavorVec GaussianTmdSet::xh1Tperp(Real x, Real Q_sq, Real k_perp_sq) const {
	return tmd_gaussian_factor(vars.h1Tperp, k_perp_sq)*xh1Tperp(x, Q_sq);
}
FlavorVec GaussianTmdSet::xh(Real x, Real Q_sq, Real k_perp_sq) const {
	return tmd_gaussian_factor(vars.h, k_perp_sq)*xh(x, Q_sq);
}
FlavorVec GaussianTmdSet::xhL(Real x, Real Q_sq, Real k_perp_sq) const {
	return tmd_gaussian_factor(vars.hL, k_perp_sq)*xhL(x, Q_sq);
}
FlavorVec GaussianTmdSet::xhT(Real x, Real Q_sq, Real k_perp_sq) const {
	return tmd_gaussian_factor(vars.hT, k_perp_sq)*xhT(x, Q_sq);
}
FlavorVec GaussianTmdSet::xhTperp(Real x, Real Q_sq, Real k_perp_sq) const {
	return tmd_gaussian_factor(vars.hTperp, k_perp_sq)*xhTperp(x, Q_sq);
}
FlavorVec GaussianTmdSet::xe(Real x, Real Q_sq, Real k_perp_sq) const {
	return tmd_gaussian_factor(vars.e, k_perp_sq)*xe(x, Q_sq);
}
FlavorVec GaussianTmdSet::xeL(Real x, Real Q_sq, Real k_perp_sq) const {
	return tmd_gaussian_factor(vars.eL, k_perp_sq)*xeL(x, Q_sq);
}
FlavorVec GaussianTmdSet::xeT(Real x, Real Q_sq, Real k_perp_sq) const {
	return tmd_gaussian_factor(vars.eT, k_perp_sq)*xeT(x, Q_sq);
}
FlavorVec GaussianTmdSet::xeTperp(Real x, Real Q_sq, Real k_perp_sq) const {
	return tmd_gaussian_factor(vars.eTperp, k_perp_sq)*xeTperp(x, Q_sq);
}

FlavorVec GaussianTmdSet::D1(part::Hadron h, Real z, Real Q_sq, Real p_perp_sq) const {
	return tmd_gaussian_factor(vars.D1, p_perp_sq)*D1(h, z, Q_sq);
}
FlavorVec GaussianTmdSet::H1perp(part::Hadron h, Real z, Real Q_sq, Real p_perp_sq) const {
	return tmd_gaussian_factor(vars.H1perp, p_perp_sq)*H1perp(h, z, Q_sq);
}
FlavorVec GaussianTmdSet::Dperp_tilde(part::Hadron h, Real z, Real Q_sq, Real p_perp_sq) const {
	return tmd_gaussian_factor(vars.Dperp_tilde, p_perp_sq)*Dperp_tilde(h, z, Q_sq);
}
FlavorVec GaussianTmdSet::H_tilde(part::Hadron h, Real z, Real Q_sq, Real p_perp_sq) const {
	return tmd_gaussian_factor(vars.H_tilde, p_perp_sq)*H_tilde(h, z, Q_sq);
}
FlavorVec GaussianTmdSet::Gperp_tilde(part::Hadron h, Real z, Real Q_sq, Real p_perp_sq) const {
	return tmd_gaussian_factor(vars.Gperp_tilde, p_perp_sq)*Gperp_tilde(h, z, Q_sq);
}
FlavorVec GaussianTmdSet::E_tilde(part::Hadron h, Real z, Real Q_sq, Real p_perp_sq) const {
	return tmd_gaussian_factor(vars.E_tilde, p_perp_sq)*E_tilde(h, z, Q_sq);
}

// WW-type approximation.
FlavorVec WwTmdSet::xf1(Real, Real, Real) const {
	return FlavorVec(flavor_count);
}
FlavorVec WwTmdSet::xf1Tperp(Real, Real, Real) const {
	return FlavorVec(flavor_count);
}
FlavorVec WwTmdSet::xg1(Real, Real, Real) const {
	return FlavorVec(flavor_count);
}
FlavorVec WwTmdSet::xg1Tperp(Real, Real, Real) const {
	return FlavorVec(flavor_count);
}
FlavorVec WwTmdSet::xh1(Real, Real, Real) const {
	return FlavorVec(flavor_count);
}
FlavorVec WwTmdSet::xh1perp(Real, Real, Real) const {
	return FlavorVec(flavor_count);
}
FlavorVec WwTmdSet::xh1Lperp(Real, Real, Real) const {
	return FlavorVec(flavor_count);
}
FlavorVec WwTmdSet::xh1Tperp(Real, Real, Real) const {
	return FlavorVec(flavor_count);
}

FlavorVec WwTmdSet::D1(part::Hadron, Real, Real, Real) const {
	return FlavorVec(flavor_count);
}
FlavorVec WwTmdSet::H1perp(part::Hadron, Real, Real, Real) const {
	return FlavorVec(flavor_count);
}

FlavorVec WwTmdSet::xf1TperpM1(Real x, Real Q_sq, Real k_perp_sq) const {
	Real M = mass(target);
	return (k_perp_sq/(2*M*M))*xf1Tperp(x, Q_sq, k_perp_sq);
}
FlavorVec WwTmdSet::xg1TperpM1(Real x, Real Q_sq, Real k_perp_sq) const {
	Real M = mass(target);
	return (k_perp_sq/(2*M*M))*xg1Tperp(x, Q_sq, k_perp_sq);
}
FlavorVec WwTmdSet::xh1perpM1(Real x, Real Q_sq, Real k_perp_sq) const {
	Real M = mass(target);
	return (k_perp_sq/(2*M*M))*xh1perp(x, Q_sq, k_perp_sq);
}
FlavorVec WwTmdSet::xh1LperpM1(Real x, Real Q_sq, Real k_perp_sq) const {
	Real M = mass(target);
	return (k_perp_sq/(2*M*M))*xh1Lperp(x, Q_sq, k_perp_sq);
}
FlavorVec WwTmdSet::xh1TperpM1(Real x, Real Q_sq, Real k_perp_sq) const {
	Real M = mass(target);
	return (k_perp_sq/(2*M*M))*xh1Tperp(x, Q_sq, k_perp_sq);
}

FlavorVec WwTmdSet::xfT(Real x, Real Q_sq, Real k_perp_sq) const {
	return -xf1TperpM1(x, Q_sq, k_perp_sq)/x;
}
FlavorVec WwTmdSet::xfperp(Real x, Real Q_sq, Real k_perp_sq) const {
	return xf1(x, Q_sq, k_perp_sq)/x;
}
FlavorVec WwTmdSet::xfLperp(Real, Real, Real) const {
	return FlavorVec(flavor_count);
}
FlavorVec WwTmdSet::xfTperp(Real x, Real Q_sq, Real k_perp_sq) const {
	return xf1Tperp(x, Q_sq, k_perp_sq)/x;
}
FlavorVec WwTmdSet::xgT(Real x, Real Q_sq, Real k_perp_sq) const {
	return xg1TperpM1(x, Q_sq, k_perp_sq)/x;
}
FlavorVec WwTmdSet::xgperp(Real, Real, Real) const {
	return FlavorVec(flavor_count);
}
FlavorVec WwTmdSet::xgLperp(Real x, Real Q_sq, Real k_perp_sq) const {
	return xg1(x, Q_sq, k_perp_sq)/x;
}
FlavorVec WwTmdSet::xgTperp(Real x, Real Q_sq, Real k_perp_sq) const {
	return xg1Tperp(x, Q_sq, k_perp_sq)/x;
}
FlavorVec WwTmdSet::xh(Real x, Real Q_sq, Real k_perp_sq) const {
	return -2.*xh1perpM1(x, Q_sq, k_perp_sq)/x;
}
FlavorVec WwTmdSet::xhL(Real x, Real Q_sq, Real k_perp_sq) const {
	return -2.*xh1LperpM1(x, Q_sq, k_perp_sq)/x;
}
FlavorVec WwTmdSet::xhT(Real x, Real Q_sq, Real k_perp_sq) const {
	return -(xh1(x, Q_sq, k_perp_sq) + xh1TperpM1(x, Q_sq, k_perp_sq))/x;
}
FlavorVec WwTmdSet::xhTperp(Real x, Real Q_sq, Real k_perp_sq) const {
	return (xh1(x, Q_sq, k_perp_sq) - xh1TperpM1(x, Q_sq, k_perp_sq))/x;
}
FlavorVec WwTmdSet::xe(Real, Real, Real) const {
	return FlavorVec(flavor_count);
}
FlavorVec WwTmdSet::xeL(Real, Real, Real) const {
	return FlavorVec(flavor_count);
}
FlavorVec WwTmdSet::xeT(Real, Real, Real) const {
	return FlavorVec(flavor_count);
}
FlavorVec WwTmdSet::xeTperp(Real, Real, Real) const {
	return FlavorVec(flavor_count);
}

FlavorVec WwTmdSet::Dperp_tilde(part::Hadron, Real, Real, Real) const {
	return FlavorVec(flavor_count);
}
FlavorVec WwTmdSet::H_tilde(part::Hadron, Real, Real, Real) const {
	return FlavorVec(flavor_count);
}
FlavorVec WwTmdSet::Gperp_tilde(part::Hadron, Real, Real, Real) const {
	return FlavorVec(flavor_count);
}
FlavorVec WwTmdSet::E_tilde(part::Hadron, Real, Real, Real) const {
	return FlavorVec(flavor_count);
}

GaussianWwTmdVars::GaussianWwTmdVars() :
	f1(1, INF),
	f1Tperp(1, INF),
	fT(1, INF),
	g1(1, INF),
	g1Tperp(1, INF),
	gT(1, INF),
	h1(1, INF),
	h1perp(1, INF),
	h1Lperp(1, INF),
	h1Tperp(1, INF),
	h(1, INF),
	hL(1, INF),
	hT(1, INF),
	hTperp(1, INF),
	D1(1, INF),
	H1perp(1, INF) { }

// Gaussian and WW-type approximation combined.
GaussianWwTmdSet::GaussianWwTmdSet(
	unsigned flavor_count,
	FlavorVec const& charges,
	part::Nucleus target,
	GaussianWwTmdVars const& vars) :
	GaussianTmdSet(flavor_count, charges, target, GaussianTmdVars(vars)) { }

FlavorVec GaussianWwTmdSet::xf1(Real, Real) const {
	return FlavorVec(flavor_count);
}
FlavorVec GaussianWwTmdSet::xf1Tperp(Real, Real) const {
	return FlavorVec(flavor_count);
}
FlavorVec GaussianWwTmdSet::xg1(Real, Real) const {
	return FlavorVec(flavor_count);
}
FlavorVec GaussianWwTmdSet::xg1Tperp(Real, Real) const {
	return FlavorVec(flavor_count);
}
FlavorVec GaussianWwTmdSet::xh1(Real, Real) const {
	return FlavorVec(flavor_count);
}
FlavorVec GaussianWwTmdSet::xh1perp(Real, Real) const {
	return FlavorVec(flavor_count);
}
FlavorVec GaussianWwTmdSet::xh1Lperp(Real, Real) const {
	return FlavorVec(flavor_count);
}
FlavorVec GaussianWwTmdSet::xh1Tperp(Real, Real) const {
	return FlavorVec(flavor_count);
}

FlavorVec GaussianWwTmdSet::D1(part::Hadron, Real, Real) const {
	return FlavorVec(flavor_count);
}
FlavorVec GaussianWwTmdSet::H1perp(part::Hadron, Real, Real) const {
	return FlavorVec(flavor_count);
}

// Following [2], we take the approach of first integrating over `k_perp_sq`
// with the Gaussian approximation, and then using Ww-type approximation. This
// is done to avoid contradictions between these two approximations (see section
// 4.4 of [2]).
FlavorVec GaussianWwTmdSet::xf1TperpM1(Real x, Real Q_sq) const {
	Real M = mass(target);
	return vars.f1Tperp*xf1Tperp(x, Q_sq)/(2*M*M);
}
FlavorVec GaussianWwTmdSet::xg1TperpM1(Real x, Real Q_sq) const {
	Real M = mass(target);
	return vars.g1Tperp*xg1Tperp(x, Q_sq)/(2.*M*M);
}
FlavorVec GaussianWwTmdSet::xh1perpM1(Real x, Real Q_sq) const {
	Real M = mass(target);
	return vars.h1perp*xh1perp(x, Q_sq)/(2.*M*M);
}
FlavorVec GaussianWwTmdSet::xh1LperpM1(Real x, Real Q_sq) const {
	Real M = mass(target);
	return vars.h1Lperp*xh1Lperp(x, Q_sq)/(2.*M*M);
}
FlavorVec GaussianWwTmdSet::xh1TperpM1(Real x, Real Q_sq) const {
	Real M = mass(target);
	return vars.h1Tperp*xh1Tperp(x, Q_sq)/(2.*M*M);
}

FlavorVec GaussianWwTmdSet::xfT(Real x, Real Q_sq) const {
	return -xf1TperpM1(x, Q_sq)/x;
}
FlavorVec GaussianWwTmdSet::xfperp(Real x, Real Q_sq) const {
	return xf1(x, Q_sq)/x;
}
FlavorVec GaussianWwTmdSet::xfLperp(Real, Real) const {
	return FlavorVec(flavor_count);
}
FlavorVec GaussianWwTmdSet::xfTperp(Real x, Real Q_sq) const {
	return xf1Tperp(x, Q_sq)/x;
}
FlavorVec GaussianWwTmdSet::xgT(Real x, Real Q_sq) const {
	return xg1TperpM1(x, Q_sq)/x;
}
FlavorVec GaussianWwTmdSet::xgperp(Real, Real) const {
	return FlavorVec(flavor_count);
}
FlavorVec GaussianWwTmdSet::xgLperp(Real x, Real Q_sq) const {
	return xg1(x, Q_sq)/x;
}
FlavorVec GaussianWwTmdSet::xgTperp(Real x, Real Q_sq) const {
	return xg1Tperp(x, Q_sq)/x;
}
FlavorVec GaussianWwTmdSet::xh(Real x, Real Q_sq) const {
	return -2.*xh1perpM1(x, Q_sq)/x;
}
FlavorVec GaussianWwTmdSet::xhL(Real x, Real Q_sq) const {
	return -2.*xh1LperpM1(x, Q_sq)/x;
}
FlavorVec GaussianWwTmdSet::xhT(Real x, Real Q_sq) const {
	return -(xh1(x, Q_sq) + xh1TperpM1(x, Q_sq))/x;
}
FlavorVec GaussianWwTmdSet::xhTperp(Real x, Real Q_sq) const {
	return (xh1(x, Q_sq) - xh1TperpM1(x, Q_sq))/x;
}
FlavorVec GaussianWwTmdSet::xe(Real, Real) const {
	return FlavorVec(flavor_count);
}
FlavorVec GaussianWwTmdSet::xeL(Real, Real) const {
	return FlavorVec(flavor_count);
}
FlavorVec GaussianWwTmdSet::xeT(Real, Real) const {
	return FlavorVec(flavor_count);
}
FlavorVec GaussianWwTmdSet::xeTperp(Real, Real) const {
	return FlavorVec(flavor_count);
}

FlavorVec GaussianWwTmdSet::Dperp_tilde(part::Hadron, Real, Real) const {
	return FlavorVec(flavor_count);
}
FlavorVec GaussianWwTmdSet::H_tilde(part::Hadron, Real, Real) const {
	return FlavorVec(flavor_count);
}
FlavorVec GaussianWwTmdSet::Gperp_tilde(part::Hadron, Real, Real) const {
	return FlavorVec(flavor_count);
}
FlavorVec GaussianWwTmdSet::E_tilde(part::Hadron, Real, Real) const {
	return FlavorVec(flavor_count);
}

