#include "sidis/sf_set/stf.hpp"

#include <algorithm>
#include <cmath>
#include <fstream>
#include <limits>
#include <string>
#include <utility>

#include <LHAPDF/LHAPDF.h>

#include "sidis/constant.hpp"
#include "sidis/extra/exception.hpp"
#include "sidis/extra/interpolate.hpp"
#include "sidis/extra/math.hpp"

#define SF_SET_DIR "sidis/sf_set"
#define STF_DIR "stf"

using namespace sidis;
using namespace sidis::interp;
using namespace sidis::math;
using namespace sidis::sf;
using namespace sidis::sf::set;

namespace {

// The following parameters come from appendix A in [2].

// 5 lightest quarks/anti-quarks, neglecting gluons.
unsigned const NUM_FLAVORS = 6;

// We only have data files for proton structure function.
Real const M = MASS_P;

const Real TMD_ZERO_WIDTH[NUM_FLAVORS] = {
	INF, INF,
	INF, INF,
	INF, INF,
};

// General parameters.
Real const TMD_SHAPE_LAMBDA = 0.04;
Real const TMD_SHAPE_Q_SQ_0 = 2.;

// Parameters for PDF.
const Real F1_MEAN_K_PERP_SQ[NUM_FLAVORS] = {
	0.50516883, 0.52059222,
	0.50516883, 0.52059222,
	0.52059222, 0.52059222,
};
const Real D1_MEAN_P_PERP_SQ[NUM_FLAVORS] = {
	0.12605078, 0.14277642,
	0.14277642, 0.12605078,
	0.14277642, 0.14277642,
};

// Parameters for transversity.
Real const H1_ALPHA = 1.11;
Real const H1_BETA = 3.64;
Real const H1_N[6] = {
	0.46, -1.00, 0.,
	0., 0., 0.,
};
const Real H1_MEAN_K_PERP_SQ[NUM_FLAVORS] = {
	0.80400386, 0.835004,
	0.80400386, 0.835004,
	0.835004, 0.835004,
};
const Real H1_SHAPE_1[NUM_FLAVORS][10] = {
	{ 0.1, 0.57058867, 4.60532016, 0., 0., -0.015, 0., 0.5, 0., 0. },
	{ 0., 1., 10., 0., 0., 0., 0., 0., 0., 0. },
	{ -0.2, 4., 4.60532016, 0., 0., 0.015, 0., 0.5, 0., 0. },
	{ 0., 1., 10., 0., 0., 0., 0., 0., 0., 0. },
	{ 0., 1., 10., 0., 0., 0., 0., 0., 0., 0. },
	{ 0., 1., 10., 0., 0., 0., 0., 0., 0., 0. },
};
const Real H1_SHAPE_2[NUM_FLAVORS][10] = { 0. };

// Parameters for Collins.
const Real COLLINS_MEAN_P_PERP_SQ[NUM_FLAVORS] = {
	0.12337425, 0.06056566,
	0.06056566, 0.12337425,
	0.06056566, 0.06056566,
};
const Real COLLINS_SHAPE_1[NUM_FLAVORS][10] = {
	{ 1.19421812e-01, 3.12545022e+00, 2.94187554e+00, 0.00000000e+00, 0.00000000e+00, -1.50000000e-03, 2.50000000e-01, 1.00000000e+00, 0.00000000e+00, 0.00000000e+00 },
	{ -2.93124180e-01, 5.00000000e-01, 2.14753393e+00, 0.00000000e+00, 0.00000000e+00, 1.50000000e-03, 1.41311020e-15, 7.54615240e-02, 0.00000000e+00, 0.00000000e+00 },
	{ -2.93124180e-01, 5.00000000e-01, 2.14753393e+00, 0.00000000e+00, 0.00000000e+00, 1.50000000e-03, 1.41311020e-15, 7.54615240e-02, 0.00000000e+00, 0.00000000e+00 },
	{ 1.19421812e-01, 3.12545022e+00, 2.94187554e+00, 0.00000000e+00, 0.00000000e+00, -1.50000000e-03, 2.50000000e-01, 1.00000000e+00, 0.00000000e+00, 0.00000000e+00 },
	{ -2.93124180e-01, 5.00000000e-01, 2.14753393e+00, 0.00000000e+00, 0.00000000e+00, 1.50000000e-03, 1.41311020e-15, 7.54615240e-02, 0.00000000e+00, 0.00000000e+00 },
	{ -2.93124180e-01, 5.00000000e-01, 2.14753393e+00, 0.00000000e+00, 0.00000000e+00, 1.50000000e-03, 1.41311020e-15, 7.54615240e-02, 0.00000000e+00, 0.00000000e+00 },
};
const Real COLLINS_SHAPE_2[NUM_FLAVORS][10] = {
	{ 1000., 0., 10., 0., 0., 0., 0., 4.24477475e-01, 0., 0. },
	{ 200., 0., 12.80882825, 0., 0., 0., 0., 0.49999953, 0., 0. },
	{ 200., 0., 12.80882825, 0., 0., 0., 0., 0.49999953, 0., 0. },
	{ 1000., 0., 10., 0., 0., 0., 0., 4.24477475e-01, 0., 0. },
	{ 200., 0., 12.80882825, 0., 0., 0., 0., 0.49999953, 0., 0. },
	{ 200., 0., 12.80882825, 0., 0., 0., 0., 0.49999953, 0., 0. },
};

// Parameters for Sivers.
const Real SIVERS_MEAN_K_PERP_SQ[NUM_FLAVORS] = {
	0.41260292, 0.488407,
	0.41260292, 0.488407,
	0.488407,   0.488407,
};
const Real SIVERS_SHAPE_1[NUM_FLAVORS][10] = {
	{ -9.18750922e-03, 5.74141750e-01, 7.44642855e+00, 0., 0., 5.00000000e-03, 0., 5.00111886e-01, 0., 0. },
	{ 0., 0., 0., 0., 0., 0., 0., 0., 0., 0. },
	{ 1.12128225e-02, 3.66897990e-01, 7.44642855e+00, 0., 0., -5.00000000e-03, 0., 5.00111886e-01, 0., 0. },
	{ 0., 0., 0., 0., 0., 0., 0., 0., 0., 0. },
	{ 0., 0., 0., 0., 0., 0., 0., 0., 0., 0. },
};
const Real SIVERS_SHAPE_2[NUM_FLAVORS][10] = { 0. };

Real tmd_shape_x(Real x, Real s, Real const* params) {
	Real N = params[0] + s*params[5];
	Real a = params[1] + s*params[6];
	Real b = params[2] + s*params[7];
	Real c = params[3] + s*params[8];
	Real d = params[4] + s*params[9];
    return N*std::pow(x, a)*std::pow(1. - x, b)*(1. + c*x + d*x*x);
	
}

Real beta(Real x, Real y) {
	return std::tgamma(x)*std::tgamma(y)/std::tgamma(x + y);
}

Real tmd_shape(Real x, Real Q_sq, Real const* params1, Real const* params2) {
	Real s = std::log(
		std::log(Q_sq/TMD_SHAPE_LAMBDA)/std::log(TMD_SHAPE_Q_SQ_0/TMD_SHAPE_LAMBDA));
	Real shape1 = tmd_shape_x(x, s, params1);
	Real shape2 = tmd_shape_x(x, s, params2);
	// TODO: Find why this condition is here.
	if (params2[0] != 0.) {
		shape1 *= 1. + shape2;
	}
	Real N1 = beta(params1[1] + 2., params1[2] + 1.)
		+ params1[3]*beta(params1[1] + 3., params1[2] + 1.)
		+ params1[4]*beta(params1[1] + 4., params1[2] + 1.);
	Real N2 = params2[0]*(
		beta(params1[1] + params2[1] + 2., params1[2] + params2[2] + 1.)
		+ params1[4]*params2[4]
			*beta(params1[1] + params2[1] + 6., params1[2] + params2[2] + 1.));
	Real N3 = params2[0]*(params1[3] + params2[3])
		*beta(params1[1] + params2[1] + 3., params1[2] + params2[2] + 1.);
	Real N4 = params2[0]*(params1[3]*params2[4] + params2[3]*params1[4])
		*beta(params1[1] + params2[1] + 5., params1[2] + params2[2] + 1.);
	Real N5 = params2[0]*(params1[3]*params2[3] + params1[4] + params2[4])
		*beta(params1[1] + params2[1] + 4., params1[2] + params2[2] + 1);
	Real N = N1 + N2 + N3 + N4 + N5;
    return shape1 / N;
}

int pdg_id(unsigned fl) {
	switch (fl) {
	case 0:
		return 2;
	case 1:
		return -2;
	case 2:
		return 1;
	case 3:
		return -1;
	case 4:
		return 3;
	case 5:
		return -3;
	default:
		return std::numeric_limits<int>::max();
	}
}

}

struct StfTmdSet::Impl {
	// PDF are calculated with the LHAPDF library.
	LHAPDF::PDF* tmd_set;
	LHAPDF::PDF* ff_set;

	Impl() {
		LHAPDF::setVerbosity(0);
		LHAPDF::pathsPrepend(DATADIR "/" SF_SET_DIR "/" STF_DIR "/");
		LHAPDF::pathsPrepend("../share/" SF_SET_DIR "/" STF_DIR "/");
		LHAPDF::pathsPrepend(SF_SET_DIR "/" STF_DIR "/");
		LHAPDF::pathsPrepend(STF_DIR "/");
		LHAPDF::pathsPrepend("./");
		tmd_set = LHAPDF::mkPDF("CJ15lo", 0);
		ff_set = LHAPDF::mkPDF("dsspipLO", 0);
		ff_set->setForcePositive(0);
	}
	Impl(Impl const& impl) = delete;
	Impl& operator=(Impl const& impl) = delete;
	~Impl() {
		delete tmd_set;
		delete ff_set;
	}
};

StfTmdSet::StfTmdSet() :
	GaussianTmdSet(
		NUM_FLAVORS,
		part::Nucleus::P,
		// `mean_f1`.
		F1_MEAN_K_PERP_SQ,
		// `mean_f1Tperp`.
		SIVERS_MEAN_K_PERP_SQ,
		// `mean_fT`.
		TMD_ZERO_WIDTH,
		// `mean_fperp`.
		TMD_ZERO_WIDTH,
		// `mean_fLperp`.
		TMD_ZERO_WIDTH,
		// `mean_fTperp`.
		TMD_ZERO_WIDTH,
		// `mean_g1`.
		TMD_ZERO_WIDTH,
		// `mean_g1Tperp`.
		TMD_ZERO_WIDTH,
		// `mean_gT`.
		TMD_ZERO_WIDTH,
		// `mean_gperp`.
		TMD_ZERO_WIDTH,
		// `mean_gLperp`.
		TMD_ZERO_WIDTH,
		// `mean_gTperp`.
		TMD_ZERO_WIDTH,
		// `mean_h1`.
		H1_MEAN_K_PERP_SQ,
		// `mean_h1perp`.
		TMD_ZERO_WIDTH,
		// `mean_h1Lperp`.
		TMD_ZERO_WIDTH,
		// `mean_h1Tperp`.
		TMD_ZERO_WIDTH,
		// `mean_h`.
		TMD_ZERO_WIDTH,
		// `mean_hL`.
		TMD_ZERO_WIDTH,
		// `mean_hT`.
		TMD_ZERO_WIDTH,
		// `mean_hTperp`.
		TMD_ZERO_WIDTH,
		// `mean_e`.
		TMD_ZERO_WIDTH,
		// `mean_eL`.
		TMD_ZERO_WIDTH,
		// `mean_eT`.
		TMD_ZERO_WIDTH,
		// `mean_eTperp`.
		TMD_ZERO_WIDTH,
		// `mean_D1`.
		D1_MEAN_P_PERP_SQ,
		// `mean_H1perp`.
		COLLINS_MEAN_P_PERP_SQ,
		// `mean_Dperp_tilde`.
		TMD_ZERO_WIDTH,
		// `mean_H_tilde`.
		TMD_ZERO_WIDTH,
		// `mean_Gperp_tilde`.
		TMD_ZERO_WIDTH,
		// `mean_E_tilde`.
		TMD_ZERO_WIDTH),
	_impl(new Impl()) { }

StfTmdSet::StfTmdSet(StfTmdSet&& other) noexcept :
		GaussianTmdSet(
			NUM_FLAVORS,
			part::Nucleus::P,
			F1_MEAN_K_PERP_SQ,
			SIVERS_MEAN_K_PERP_SQ,
			TMD_ZERO_WIDTH,
			TMD_ZERO_WIDTH,
			TMD_ZERO_WIDTH,
			TMD_ZERO_WIDTH,
			TMD_ZERO_WIDTH,
			TMD_ZERO_WIDTH,
			TMD_ZERO_WIDTH,
			TMD_ZERO_WIDTH,
			TMD_ZERO_WIDTH,
			TMD_ZERO_WIDTH,
			H1_MEAN_K_PERP_SQ,
			TMD_ZERO_WIDTH,
			TMD_ZERO_WIDTH,
			TMD_ZERO_WIDTH,
			TMD_ZERO_WIDTH,
			TMD_ZERO_WIDTH,
			TMD_ZERO_WIDTH,
			TMD_ZERO_WIDTH,
			TMD_ZERO_WIDTH,
			TMD_ZERO_WIDTH,
			TMD_ZERO_WIDTH,
			TMD_ZERO_WIDTH,
			D1_MEAN_P_PERP_SQ,
			COLLINS_MEAN_P_PERP_SQ,
			TMD_ZERO_WIDTH,
			TMD_ZERO_WIDTH,
			TMD_ZERO_WIDTH,
			TMD_ZERO_WIDTH),
		_impl(nullptr) {
	std::swap(_impl, other._impl);
}
StfTmdSet& StfTmdSet::operator=(StfTmdSet&& other) noexcept {
	std::swap(_impl, other._impl);
	return *this;
}
StfTmdSet::~StfTmdSet() {
	if (_impl != nullptr) {
		delete _impl;
	}
}

Real StfTmdSet::charge(unsigned fl) const {
	switch (fl) {
	// Up.
	case 0:
		return 2./3.;
	// Up bar.
	case 1:
		return -2./3.;
	// Down.
	case 2:
		return -1./3.;
	// Down bar.
	case 3:
		return 1./3.;
	// Strange.
	case 4:
		return -1./3.;
	// Strange bar.
	case 5:
		return 1./3.;
	default:
		return 0.;
	}
}

Real StfTmdSet::xf1(unsigned fl, Real x, Real Q_sq) const {
	int pdg = pdg_id(fl);
	return _impl->tmd_set->xfxQ2(pdg, x, Q_sq);
}

Real StfTmdSet::xf1Tperp(unsigned fl, Real x, Real Q_sq) const {
	return 2.*x*M*M/SIVERS_MEAN_K_PERP_SQ[fl]
		*tmd_shape(x, Q_sq, SIVERS_SHAPE_1[fl], SIVERS_SHAPE_2[fl]);
}

Real StfTmdSet::xh1(unsigned fl, Real x, Real Q_sq) const {
	return x*tmd_shape(x, Q_sq, H1_SHAPE_1[fl], H1_SHAPE_2[fl]);
}

Real StfTmdSet::D1(part::Hadron h, unsigned fl, Real z, Real Q_sq) const {
	int pdg = pdg_id(fl);
	switch (h) {
	case part::Hadron::PI_P:
		return _impl->ff_set->xfxQ2(pdg, z, Q_sq)/z;
	case part::Hadron::PI_M:
		// TODO: Handle this case later.
	default:
		throw HadronOutOfRange(h);
	}
}

Real StfTmdSet::H1perp(part::Hadron h, unsigned fl, Real z, Real Q_sq) const {
	Real mh = mass(h);
	switch (h) {
	case part::Hadron::PI_P:
		return 2.*z*z*mh*mh/COLLINS_MEAN_P_PERP_SQ[fl]
			*tmd_shape(z, Q_sq, COLLINS_SHAPE_1[fl], COLLINS_SHAPE_2[fl]);
	case part::Hadron::PI_M:
		// TODO: Handle this case later.
	default:
		throw HadronOutOfRange(h);
	}
}

