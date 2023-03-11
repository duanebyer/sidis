#include "sidis/sf_set/prokudin.hpp"

#include <array>
#include <cmath>
#include <fstream>
#include <sstream>
#include <string>
#include <utility>
#include <vector>

#include <mstwpdf.h>

#include "sidis/constant.hpp"
#include "sidis/extra/exception.hpp"
#include "sidis/extra/interpolate.hpp"
#include "sidis/extra/math.hpp"

#define SF_SET_DIR "sidis/sf_set"
#define PROKUDIN_DIR "prokudin"

using namespace sidis;
using namespace sidis::interp;
using namespace sidis::math;
using namespace sidis::sf;
using namespace sidis::sf::set;

namespace {

// The following parameters come from appendix A in [2].

unsigned const NUM_FLAVORS = 6;

// We only have data files for proton structure function.
Real const M = MASS_P;

// Parameters for PDF.
Real const F1_MEAN_K_PERP_SQ = 0.25;
Real const G1_MEAN_K_PERP_SQ = 0.76 * F1_MEAN_K_PERP_SQ;
Real const D1_MEAN_P_PERP_SQ = 0.2;

// Parameters for transversity.
Real const H1_ALPHA = 1.11;
Real const H1_BETA = 3.64;
FlavorVec const H1_N = {
	0.46, -1.00, 0.,
	0., 0., 0.,
};
Real const H1_MEAN_K_PERP_SQ = 0.25;
// Reference [2] is not consistent about the choice of `mean_hT`, so here we
// choose something in between the two options.
Real const HT_MEAN_K_PERP_SQ = 0.15;

// Parameters for Collins.
Real const COLLINS_M_SQ = 1.50;
Real const COLLINS_M = std::sqrt(COLLINS_M_SQ);
Real const COLLINS_GAMMA = 1.06;
Real const COLLINS_DELTA = 0.07;
Real const COLLINS_N_FAV = 0.49;
Real const COLLINS_N_DISFAV = -1.;
Real const COLLINS_MEAN_P_PERP_SQ = (D1_MEAN_P_PERP_SQ * COLLINS_M_SQ)
	/ (D1_MEAN_P_PERP_SQ + COLLINS_M_SQ);

// Parameters for Sivers.
Real const SIVERS_M_1_SQ = 0.19;
Real const SIVERS_M_1 = std::sqrt(SIVERS_M_1_SQ);
FlavorVec const SIVERS_ALPHA = {
	0.35, 0.44, 0.,
	0., 0., 0.,
};
FlavorVec const SIVERS_BETA = {
	2.6, 0.90, 0.,
	0., 0., 0.,
};
FlavorVec const SIVERS_N = {
	0.40, -0.97, 0.,
	0., 0., 0.,
};
Real const SIVERS_MEAN_K_PERP_SQ = (F1_MEAN_K_PERP_SQ * SIVERS_M_1_SQ)
	/ (F1_MEAN_K_PERP_SQ + SIVERS_M_1_SQ);

// Parameters for Boer-Mulders.
Real const BM_M_1_SQ = 0.34;
Real const BM_M_1 = std::sqrt(BM_M_1_SQ);
FlavorVec const BM_ALPHA = {
	0.73, 1.08, 0.79,
	0.79, 0.79, 0.79,
};
FlavorVec const BM_BETA = {
	3.46, 3.46, 3.46,
	3.46, 3.46, 3.46,
};
FlavorVec const BM_A = {
	0.35, -0.90, 0.24,
	0.04, -0.40, -1.,
};
FlavorVec const BM_LAMBDA = {
	2.1, -1.111, 0.,
	0., 0., 0.,
};
Real const BM_MEAN_K_PERP_SQ = (F1_MEAN_K_PERP_SQ * BM_M_1_SQ)
	/ (F1_MEAN_K_PERP_SQ + BM_M_1_SQ);

// Parameters for pretzelosity.
Real const PRETZ_M_TT_SQ = 0.18;
Real const PRETZ_M_TT = std::sqrt(PRETZ_M_TT_SQ);
Real const PRETZ_ALPHA = 2.5;
Real const PRETZ_BETA = 2.;
FlavorVec const PRETZ_N = {
	1., -1., 0.,
	0., 0., 0.,
};
Real const PRETZ_MEAN_K_PERP_SQ = (F1_MEAN_K_PERP_SQ * PRETZ_M_TT_SQ)
	/ (F1_MEAN_K_PERP_SQ + PRETZ_M_TT_SQ);

Real G(Real ph_t_sq, Real l) {
	// Equation [2.5.2].
	return std::exp(-ph_t_sq / l) / (PI * l);
}
Real lambda(Real z, Real mean_kperp_sq, Real mean_pperp_sq) {
	return sq(z) * mean_kperp_sq + mean_pperp_sq;
}

// Struct summarizing the TMD Gaussian variances.
GaussianWwTmdVars const TMD_VARS = []() {
	GaussianWwTmdVars vars{};
	vars.f1 = FlavorVec(F1_MEAN_K_PERP_SQ);
	vars.f1Tperp = FlavorVec(SIVERS_MEAN_K_PERP_SQ);
	// vars.fT
	vars.g1 = FlavorVec(G1_MEAN_K_PERP_SQ);
	vars.g1Tperp = FlavorVec(G1_MEAN_K_PERP_SQ);
	vars.gT = FlavorVec(G1_MEAN_K_PERP_SQ);
	vars.h1 = FlavorVec(H1_MEAN_K_PERP_SQ);
	vars.h1perp = FlavorVec(BM_MEAN_K_PERP_SQ);
	vars.h1Lperp = FlavorVec(H1_MEAN_K_PERP_SQ);
	vars.h1Tperp = FlavorVec(PRETZ_MEAN_K_PERP_SQ);
	// vars.h
	vars.hL = FlavorVec(H1_MEAN_K_PERP_SQ);
	vars.hT = FlavorVec(HT_MEAN_K_PERP_SQ);
	vars.hTperp = FlavorVec(HT_MEAN_K_PERP_SQ);
	vars.D1 = FlavorVec(D1_MEAN_P_PERP_SQ);
	vars.H1perp = FlavorVec(COLLINS_MEAN_P_PERP_SQ);
	// `fT` and `h` are neglected because of a lack of information about the
	// `fT` TMD (see the discussion in section 7.6 of [2]).
	return vars;
}();

FlavorVec const CHARGES = {
	+2./3., // Up.
	-1./3., // Down.
	-1./3., // Strange.
	-2./3., // Up bar.
	+1./3., // Down bar.
	+1./3., // Strange bar.
};

// Finds a grid file.
std::istream& find_file(std::ifstream& fin, char const* file_name) {
	fin.open(
		std::string(DATADIR "/" SF_SET_DIR "/" PROKUDIN_DIR "/") + file_name);
	if (fin) {
		return fin;
	}
	fin.open(std::string("../share/" SF_SET_DIR "/" PROKUDIN_DIR "/") + file_name);
	if (fin) {
		return fin;
	}
	fin.open(std::string(SF_SET_DIR "/" PROKUDIN_DIR "/") + file_name);
	if (fin) {
		return fin;
	}
	fin.open(std::string(PROKUDIN_DIR "/") + file_name);
	if (fin) {
		return fin;
	}
	fin.open(file_name);
	if (fin) {
		return fin;
	}
	throw DataFileNotFound(file_name);
}

// Convenience method for loading grid data from a file. Automatically searches
// in several possible directories for the file.
template<typename T, std::size_t N, std::size_t K>
std::array<Grid<T, N>, K> load_grids(char const* file_name) {
	std::ifstream in;
	find_file(in, file_name);
	std::vector<std::array<T, N + K> > data;
	std::string line;
	while (std::getline(in, line)) {
		std::stringstream line_in(line);
		std::array<T, N + K> next;
		for (std::size_t idx = 0; idx < N + K; ++idx) {
			line_in >> next[idx];
		}
		data.push_back(next);
		if (!line_in) {
			throw DataFileParseError(file_name);
		}
	}
	// By default, assume the grids are only accurate to single precision.
	return read_grids<T, N, K>(data, 0.000001);
}

// Shared implementation between `ProkudinTmdSet` and `ProkudinSfSet`.
struct ProkudinImpl {
	// PDF are calculated with the MSTWPDF library.
	std::ifstream file_pdf;
	mstw::c_mstwpdf pdf;

	// Load data files from WWSIDIS repository for the TMDs and FFs.
	std::array<Grid<Real, 2>, 6> data_D1_pi_plus;
	std::array<Grid<Real, 2>, 6> data_D1_pi_minus;
	std::array<Grid<Real, 2>, 6> data_g1;
	std::array<Grid<Real, 2>, 6> data_xgT;
	std::array<Grid<Real, 2>, 2> data_xh1LperpM1;
	// Soffer bound.
	std::array<Grid<Real, 2>, 6> data_sb;

	// Fragmentation functions.
	CubicView<Real, 2> interp_D1_pi_plus[6];
	CubicView<Real, 2> interp_D1_pi_minus[6];

	// Transverse momentum distributions.
	CubicView<Real, 2> interp_g1[6];
	CubicView<Real, 2> interp_xgT[6];
	CubicView<Real, 2> interp_xh1LperpM1[2];
	CubicView<Real, 2> interp_sb[6];

	ProkudinImpl() :
			file_pdf(),
			pdf(find_file(file_pdf, "mstw2008lo.00.dat"), false, true),
			data_D1_pi_plus(load_grids<Real, 2, 6>("fragmentationpiplus.dat")),
			data_D1_pi_minus(load_grids<Real, 2, 6>("fragmentationpiminus.dat")),
			data_g1(load_grids<Real, 2, 6>("g1.dat")),
			data_xgT(load_grids<Real, 2, 6>("gT_u_d_ubar_dbar_s_sbar.dat")),
			data_xh1LperpM1(load_grids<Real, 2, 2>("xh1Lperp_u_d.dat")),
			data_sb(load_grids<Real, 2, 6>("SofferBound.dat")),
			interp_D1_pi_plus{
				CubicView<Real, 2>(data_D1_pi_plus[0]),
				CubicView<Real, 2>(data_D1_pi_plus[1]),
				CubicView<Real, 2>(data_D1_pi_plus[2]),
				CubicView<Real, 2>(data_D1_pi_plus[3]),
				CubicView<Real, 2>(data_D1_pi_plus[4]),
				CubicView<Real, 2>(data_D1_pi_plus[5]),
			},
			interp_D1_pi_minus{
				CubicView<Real, 2>(data_D1_pi_minus[0]),
				CubicView<Real, 2>(data_D1_pi_minus[1]),
				CubicView<Real, 2>(data_D1_pi_minus[2]),
				CubicView<Real, 2>(data_D1_pi_minus[3]),
				CubicView<Real, 2>(data_D1_pi_minus[4]),
				CubicView<Real, 2>(data_D1_pi_minus[5]),
			},
			interp_g1{
				CubicView<Real, 2>(data_g1[0]),
				CubicView<Real, 2>(data_g1[1]),
				CubicView<Real, 2>(data_g1[2]),
				CubicView<Real, 2>(data_g1[3]),
				CubicView<Real, 2>(data_g1[4]),
				CubicView<Real, 2>(data_g1[5]),
			},
			interp_xgT{
				CubicView<Real, 2>(data_xgT[0]),
				CubicView<Real, 2>(data_xgT[1]),
				CubicView<Real, 2>(data_xgT[4]),
				CubicView<Real, 2>(data_xgT[2]),
				CubicView<Real, 2>(data_xgT[3]),
				CubicView<Real, 2>(data_xgT[5]),
			},
			interp_xh1LperpM1{
				CubicView<Real, 2>(data_xh1LperpM1[0]),
				CubicView<Real, 2>(data_xh1LperpM1[1]),
			},
			interp_sb{
				CubicView<Real, 2>(data_sb[0]),
				CubicView<Real, 2>(data_sb[1]),
				CubicView<Real, 2>(data_sb[2]),
				CubicView<Real, 2>(data_sb[3]),
				CubicView<Real, 2>(data_sb[4]),
				CubicView<Real, 2>(data_sb[5]),
			} {
		file_pdf.close();
	}
};

FlavorVec eval_interp_arr(const CubicView<Real, 2> (&interp)[6], Real x, Real Q_sq) {
	return FlavorVec {
		interp[0]({ (x), (Q_sq) }),
		interp[1]({ (x), (Q_sq) }),
		interp[2]({ (x), (Q_sq) }),
		interp[3]({ (x), (Q_sq) }),
		interp[4]({ (x), (Q_sq) }),
		interp[5]({ (x), (Q_sq) }),
	};
}

}

// Can't use typedefs here, so use inheritance to make the implementations the
// same for both `ProkudinTmdSet` and `ProkudinSfSet`.
struct ProkudinTmdSet::Impl : public ProkudinImpl { };
struct ProkudinSfSet::Impl : public ProkudinImpl { };

ProkudinTmdSet::ProkudinTmdSet() :
	GaussianWwTmdSet(part::Nucleus::P, NUM_FLAVORS, CHARGES, TMD_VARS),
	_impl(pimpl::make_pimpl<Impl>()) { }

FlavorVec ProkudinTmdSet::xf1(Real x, Real Q_sq) const {
	Real Q = std::sqrt(Q_sq);
	return FlavorVec {
		_impl->pdf.parton(8, x, Q) + _impl->pdf.parton(-2, x, Q),
		_impl->pdf.parton(7, x, Q) + _impl->pdf.parton(-1, x, Q),
		_impl->pdf.parton(3, x, Q),
		_impl->pdf.parton(-2, x, Q),
		_impl->pdf.parton(-1, x, Q),
		_impl->pdf.parton(-3, x, Q),
	};
}

FlavorVec ProkudinTmdSet::xf1Tperp(Real x, Real Q_sq) const {
	// Equation [2.A.2].
	return -std::sqrt(2.*E)*M/SIVERS_M_1
		*SIVERS_MEAN_K_PERP_SQ/F1_MEAN_K_PERP_SQ
		*SIVERS_N
		*pow_vec(x, SIVERS_ALPHA)*pow_vec(1. - x, SIVERS_BETA)
		*pow_vec(SIVERS_ALPHA + SIVERS_BETA, SIVERS_ALPHA + SIVERS_BETA)
		*pow_vec(SIVERS_ALPHA, -SIVERS_ALPHA)
		*pow_vec(SIVERS_BETA, -SIVERS_BETA)
		*xf1(x, Q_sq);
}

FlavorVec ProkudinTmdSet::xg1(Real x, Real Q_sq) const {
	return x*eval_interp_arr(_impl->interp_g1, x, Q_sq);
}

FlavorVec ProkudinTmdSet::xg1Tperp(Real x, Real Q_sq) const {
	// We only have a grid for `gT`, so use the reverse WW-type approximation to
	// get `g1Tperp`.
	return (2.*M*M/G1_MEAN_K_PERP_SQ)*x
		*eval_interp_arr(_impl->interp_xgT, x, Q_sq);
}

FlavorVec ProkudinTmdSet::xh1(Real x, Real Q_sq) const {
	// Use the Soffer bound to get an upper limit on transversity (Equation
	// [2.A.7]).
	return x*H1_N
		*std::pow(x, H1_ALPHA)*std::pow(1. - x, H1_BETA)
		*std::pow(H1_ALPHA + H1_BETA, H1_ALPHA + H1_BETA)
		*std::pow(H1_ALPHA, -H1_ALPHA)
		*std::pow(H1_BETA, -H1_BETA)
		*eval_interp_arr(_impl->interp_sb, x, Q_sq);
}

FlavorVec ProkudinTmdSet::xh1perp(Real x, Real Q_sq) const {
	// Equation [2.A.18].
	return -std::sqrt(2.*E)*M/BM_M_1
		*BM_MEAN_K_PERP_SQ/F1_MEAN_K_PERP_SQ
		*BM_LAMBDA * BM_A
		*pow_vec(x, BM_ALPHA)*pow_vec(1. - x, BM_BETA)
		*pow_vec(BM_ALPHA + BM_BETA, BM_ALPHA + BM_BETA)
		*pow_vec(BM_ALPHA, -BM_ALPHA)
		*pow_vec(BM_BETA, -BM_BETA)
		*xf1(x, Q_sq);
}

FlavorVec ProkudinTmdSet::xh1Lperp(Real x, Real Q_sq) const {
	return FlavorVec {
		2.*sq(M)/H1_MEAN_K_PERP_SQ*_impl->interp_xh1LperpM1[0]({ x, Q_sq }),
		2.*sq(M)/H1_MEAN_K_PERP_SQ*_impl->interp_xh1LperpM1[1]({ x, Q_sq }),
		0., 0., 0., 0.,
	};
}

FlavorVec ProkudinTmdSet::xh1Tperp(Real x, Real Q_sq) const {
	// Equation [2.A.24].
	return E*sq(M)/PRETZ_M_TT_SQ
		*PRETZ_MEAN_K_PERP_SQ/F1_MEAN_K_PERP_SQ
		*PRETZ_N
		*std::pow(x, PRETZ_ALPHA) * std::pow(1. - x, PRETZ_BETA)
		*std::pow(PRETZ_ALPHA + PRETZ_BETA, PRETZ_ALPHA + PRETZ_BETA)
		*std::pow(PRETZ_ALPHA, -PRETZ_ALPHA)
		*std::pow(PRETZ_BETA, -PRETZ_BETA)
		*(xf1(x, Q_sq) - xg1(x, Q_sq));
}

FlavorVec ProkudinTmdSet::D1(part::Hadron h, Real z, Real Q_sq) const {
	switch (h) {
	case part::Hadron::PI_P:
		return eval_interp_arr(_impl->interp_D1_pi_plus, z, Q_sq);
	case part::Hadron::PI_M:
		return eval_interp_arr(_impl->interp_D1_pi_minus, z, Q_sq);
	default:
		throw HadronOutOfRange(h);
	}
}

FlavorVec ProkudinTmdSet::H1perp(part::Hadron h, Real z, Real Q_sq) const {
	Real mh = mass(h);
	FlavorVec collins_coeff(6, 0.);
	// Favored or dis-favored depending on charge of the quark.
	switch (h) {
	case part::Hadron::PI_P:
		// Up or anti-down favored, down or anti-up disfavored.
		collins_coeff = FlavorVec {
			COLLINS_N_FAV, COLLINS_N_DISFAV, 0.,
			COLLINS_N_DISFAV, COLLINS_N_FAV, 0.,
		};
		break;
	case part::Hadron::PI_M:
		// Opposite.
		collins_coeff = FlavorVec {
			COLLINS_N_DISFAV, COLLINS_N_FAV, 0.,
			COLLINS_N_FAV, COLLINS_N_DISFAV, 0.,
		};
		break;
	default:
		throw HadronOutOfRange(h);
	}
	return std::sqrt(2.*E)*z*mh/COLLINS_M
		*COLLINS_MEAN_P_PERP_SQ/D1_MEAN_P_PERP_SQ
		*collins_coeff
		*std::pow(z, COLLINS_GAMMA)*std::pow(1. - z, COLLINS_DELTA)
		*std::pow(COLLINS_GAMMA + COLLINS_DELTA, COLLINS_GAMMA + COLLINS_DELTA)
		*std::pow(COLLINS_GAMMA, -COLLINS_GAMMA)
		*std::pow(COLLINS_DELTA, -COLLINS_DELTA)
		*D1(h, z, Q_sq);
}

ProkudinSfSet::ProkudinSfSet() :
	SfSet(part::Nucleus::P),
	_impl(pimpl::make_pimpl<Impl>()) { }

Real ProkudinSfSet::F_UUT(part::Hadron h, Real x, Real z, Real Q_sq, Real ph_t_sq) const {
	// Equation [2.5.1a].
	Real result = (sq_vec(CHARGES)*xf1(x, Q_sq)*D1(h, z, Q_sq)).sum();
	Real l = lambda(z, F1_MEAN_K_PERP_SQ, D1_MEAN_P_PERP_SQ);
	return G(ph_t_sq, l)*result;
}
Real ProkudinSfSet::F_UU_cos_phih(part::Hadron h, Real x, Real z, Real Q_sq, Real ph_t_sq) const {
	// Equation [2.7.9a].
	Real Q = std::sqrt(Q_sq);
	Real ph_t = std::sqrt(ph_t_sq);
	// Uses a WW-type approximation to rewrite in terms of `xf1`.
	Real result = (sq_vec(CHARGES)*xf1(x, Q_sq)*D1(h, z, Q_sq)).sum();
	Real l = lambda(z, F1_MEAN_K_PERP_SQ, D1_MEAN_P_PERP_SQ);
	return -2.*F1_MEAN_K_PERP_SQ/Q*ph_t*(z/l)*G(ph_t_sq, l)*result;
}
Real ProkudinSfSet::F_UU_cos_2phih(part::Hadron h, Real x, Real z, Real Q_sq, Real ph_t_sq) const {
	// Equation [2.5.9a].
	Real mh = mass(h);
	Real result = (sq_vec(CHARGES)*xh1perpM1(x, Q_sq)*H1perpM1(h, z, Q_sq)).sum();
	Real l = lambda(z, BM_MEAN_K_PERP_SQ, COLLINS_MEAN_P_PERP_SQ);
	return 4.*M*mh*ph_t_sq*sq(z/l)*G(ph_t_sq, l)*result;
}

Real ProkudinSfSet::F_UL_sin_phih(part::Hadron h, Real x, Real z, Real Q_sq, Real ph_t_sq) const {
	// Equation [2.7.6a].
	Real mh = mass(h);
	Real Q = std::sqrt(Q_sq);
	Real ph_t = std::sqrt(ph_t_sq);
	// Use WW-type approximation to rewrite in terms of `xh1LperpM1`.
	Real result = (sq_vec(CHARGES)*xh1LperpM1(x, Q_sq)*H1perpM1(h, z, Q_sq)).sum();
	// Approximate width with `H1_MEAN_K_PERP_SQ`.
	Real l = lambda(z, H1_MEAN_K_PERP_SQ, COLLINS_MEAN_P_PERP_SQ);
	return -8.*M*mh*z*ph_t/(Q*l)*G(ph_t_sq, l)*result;
}
Real ProkudinSfSet::F_UL_sin_2phih(part::Hadron h, Real x, Real z, Real Q_sq, Real ph_t_sq) const {
	// Equation [2.6.2a].
	Real mh = mass(h);
	Real result = (sq_vec(CHARGES)*xh1LperpM1(x, Q_sq)*H1perpM1(h, z, Q_sq)).sum();
	// Approximate width with `H1_MEAN_K_PERP_SQ`.
	Real l = lambda(z, H1_MEAN_K_PERP_SQ, COLLINS_MEAN_P_PERP_SQ);
	return 4.*M*mh*ph_t_sq*sq(z/l)*G(ph_t_sq, l)*result;
}

Real ProkudinSfSet::F_UTT_sin_phih_m_phis(part::Hadron h, Real x, Real z, Real Q_sq, Real ph_t_sq) const {
	// Equation [2.5.7a].
	Real ph_t = std::sqrt(ph_t_sq);
	Real result = (sq_vec(CHARGES)*xf1TperpM1(x, Q_sq)*D1(h, z, Q_sq)).sum();
	Real l = lambda(z, SIVERS_MEAN_K_PERP_SQ, D1_MEAN_P_PERP_SQ);
	return -2.*M*z*ph_t/l*G(ph_t_sq, l)*result;
}
Real ProkudinSfSet::F_UT_sin_2phih_m_phis(part::Hadron h, Real x, Real z, Real Q_sq, Real ph_t_sq) const {
	// Equation [2.7.8a].
	Real mh = mass(h);
	Real Q = std::sqrt(Q_sq);
	// Use WW-type approximation to rewrite in terms of `xf1TperpM1`.
	Real result_1 = (sq_vec(CHARGES)*xf1TperpM1(x, Q_sq)*D1(h, z, Q_sq)).sum();
	// Use WW-type approximation to rewrite in terms of `h1TperpM2`.
	Real result_2 = (sq_vec(CHARGES)*xh1TperpM2(x, Q_sq)*H1perpM1(h, z, Q_sq)).sum();
	// Approximate width with `SIVERS_MEAN_K_PERP_SQ`.
	Real l_1 = lambda(z, SIVERS_MEAN_K_PERP_SQ, D1_MEAN_P_PERP_SQ);
	// Approximate width with `PRETZ_MEAN_K_PERP_SQ`.
	Real l_2 = lambda(z, PRETZ_MEAN_K_PERP_SQ, COLLINS_MEAN_P_PERP_SQ);
	// The paragraph following [2.7.8a] has a mistake in the WW-type
	// approximation linking `h1TM1 + h1TperpM1` with `h1TperpM2`, due to a
	// missing factor of 2.
	return 2.*M*ph_t_sq/Q*(
		SIVERS_MEAN_K_PERP_SQ*sq(z/l_1)*G(ph_t_sq, l_1)*result_1
		- 2.*M*mh*sq(z/l_2)*G(ph_t_sq, l_2)*result_2);
}
Real ProkudinSfSet::F_UT_sin_3phih_m_phis(part::Hadron h, Real x, Real z, Real Q_sq, Real ph_t_sq) const {
	// Equation [2.5.10a].
	Real mh = mass(h);
	Real ph_t = std::sqrt(ph_t_sq);
	Real result = (sq_vec(CHARGES)*xh1TperpM2(x, Q_sq)*H1perpM1(h, z, Q_sq)).sum();
	Real l = lambda(z, PRETZ_MEAN_K_PERP_SQ, COLLINS_MEAN_P_PERP_SQ);
	return 2.*sq(M)*mh*std::pow(z*ph_t/l, 3)*G(ph_t_sq, l)*result;
}
Real ProkudinSfSet::F_UT_sin_phis(part::Hadron h, Real x, Real z, Real Q_sq, Real ph_t_sq) const {
	// Equation [2.7.7a].
	Real mh = mass(h);
	Real Q = std::sqrt(Q_sq);
	// WW-type approximation used here (see [2] for details).
	Real result = (sq_vec(CHARGES)*xh1M1(x, Q_sq)*H1perpM1(h, z, Q_sq)).sum();
	Real l = lambda(z, PRETZ_MEAN_K_PERP_SQ, COLLINS_MEAN_P_PERP_SQ);
	return 8.*sq(M)*mh*sq(z)/(Q*l)*(1. - ph_t_sq/l)*G(ph_t_sq, l)*result;
}
Real ProkudinSfSet::F_UT_sin_phih_p_phis(part::Hadron h, Real x, Real z, Real Q_sq, Real ph_t_sq) const {
	// Equation [2.5.8a].
	Real mh = mass(h);
	Real ph_t = std::sqrt(ph_t_sq);
	Real result = (sq_vec(CHARGES)*xh1(x, Q_sq)*H1perpM1(h, z, Q_sq)).sum();
	Real l = lambda(z, H1_MEAN_K_PERP_SQ, COLLINS_MEAN_P_PERP_SQ);
	return 2.*mh*z*ph_t/l*G(ph_t_sq, l)*result;
}

Real ProkudinSfSet::F_LL(part::Hadron h, Real x, Real z, Real Q_sq, Real ph_t_sq) const {
	// Equation [2.5.5a].
	Real result = (sq_vec(CHARGES)*xg1(x, Q_sq)*D1(h, z, Q_sq)).sum();
	Real l = lambda(z, G1_MEAN_K_PERP_SQ, D1_MEAN_P_PERP_SQ);
	return G(ph_t_sq, l)*result;
}
Real ProkudinSfSet::F_LL_cos_phih(part::Hadron h, Real x, Real z, Real Q_sq, Real ph_t_sq) const {
	// Equation [2.7.5a].
	Real Q = std::sqrt(Q_sq);
	Real ph_t = std::sqrt(ph_t_sq);
	// Uses a WW-type approximation to rewrite in terms of `xg1`.
	Real result = (sq_vec(CHARGES)*xg1(x, Q_sq)*D1(h, z, Q_sq)).sum();
	// Approximate width with `G1_MEAN_K_PERP_SQ`.
	Real l = lambda(z, G1_MEAN_K_PERP_SQ, D1_MEAN_P_PERP_SQ);
	return -2.*G1_MEAN_K_PERP_SQ*z*ph_t/(Q*l)*G(ph_t_sq, l)*result;
}

Real ProkudinSfSet::F_LT_cos_phih_m_phis(part::Hadron h, Real x, Real z, Real Q_sq, Real ph_t_sq) const {
	// Equation [2.6.1a].
	Real ph_t = std::sqrt(ph_t_sq);
	// Uses a WW-type approximation to rewrite in terms of `xgT`.
	Real result = (sq_vec(CHARGES)*xgT(x, Q_sq)*D1(h, z, Q_sq)).sum();
	// Approximate width with `G1_MEAN_K_PERP_SQ`.
	Real l = lambda(z, G1_MEAN_K_PERP_SQ, D1_MEAN_P_PERP_SQ);
	return 2.*M*x*z*ph_t/l*G(ph_t_sq, l)*result;
}
Real ProkudinSfSet::F_LT_cos_2phih_m_phis(part::Hadron h, Real x, Real z, Real Q_sq, Real ph_t_sq) const {
	// Equation [2.7.4a].
	Real Q = std::sqrt(Q_sq);
	// Uses a WW-type approximation to rewrite in terms of `xgT`.
	Real result = (sq_vec(CHARGES)*xgT(x, Q_sq)*D1(h, z, Q_sq)).sum();
	// Approximate width with `G1_MEAN_K_PERP_SQ`.
	Real l = lambda(z, G1_MEAN_K_PERP_SQ, D1_MEAN_P_PERP_SQ);
	return -2.*G1_MEAN_K_PERP_SQ*M*x*ph_t_sq*sq(z/l)/Q*G(ph_t_sq, l)*result;
}
Real ProkudinSfSet::F_LT_cos_phis(part::Hadron h, Real x, Real z, Real Q_sq, Real ph_t_sq) const {
	// Equation [2.7.2a].
	Real Q = std::sqrt(Q_sq);
	Real result = (sq_vec(CHARGES)*xgT(x, Q_sq)*D1(h, z, Q_sq)).sum();
	// Approximate width with `G1_MEAN_K_PERP_SQ`.
	Real l = lambda(z, G1_MEAN_K_PERP_SQ, D1_MEAN_P_PERP_SQ);
	return -2.*M*x/Q*G(ph_t_sq, l)*result;
}

// Fragmentation functions.
FlavorVec ProkudinSfSet::D1(part::Hadron h, Real z, Real Q_sq) const {
	switch (h) {
	case part::Hadron::PI_P:
		return eval_interp_arr(_impl->interp_D1_pi_plus, z, Q_sq);
	case part::Hadron::PI_M:
		return eval_interp_arr(_impl->interp_D1_pi_minus, z, Q_sq);
	default:
		throw HadronOutOfRange(h);
	}
}
FlavorVec ProkudinSfSet::H1perpM1(part::Hadron h, Real z, Real Q_sq) const {
	Real mh = mass(h);
	FlavorVec collins_coeff(6, 0.);
	switch (h) {
	case part::Hadron::PI_P:
		collins_coeff = FlavorVec {
			COLLINS_N_FAV, COLLINS_N_DISFAV, 0.,
			COLLINS_N_DISFAV, COLLINS_N_FAV, 0.,
		};
		break;
	case part::Hadron::PI_M:
		collins_coeff = FlavorVec {
			COLLINS_N_DISFAV, COLLINS_N_FAV, 0.,
			COLLINS_N_FAV, COLLINS_N_DISFAV, 0.,
		};
		break;
	default:
		throw HadronOutOfRange(h);
	}
	return std::sqrt(E/2.)/(z*mh*COLLINS_M)
		*sq(COLLINS_MEAN_P_PERP_SQ)/D1_MEAN_P_PERP_SQ
		*collins_coeff
		*std::pow(z, COLLINS_GAMMA)*std::pow(1. - z, COLLINS_DELTA)
		*std::pow(COLLINS_GAMMA + COLLINS_DELTA, COLLINS_GAMMA + COLLINS_DELTA)
		*std::pow(COLLINS_GAMMA, -COLLINS_GAMMA)
		*std::pow(COLLINS_DELTA, -COLLINS_DELTA)
		*D1(h, z, Q_sq);
}

// Parton distribution functions.
FlavorVec ProkudinSfSet::xf1(Real x, Real Q_sq) const {
	Real Q = std::sqrt(Q_sq);
	return FlavorVec {
		_impl->pdf.parton(8, x, Q) + _impl->pdf.parton(-2, x, Q),
		_impl->pdf.parton(7, x, Q) + _impl->pdf.parton(-1, x, Q),
		_impl->pdf.parton(3, x, Q),
		_impl->pdf.parton(-2, x, Q),
		_impl->pdf.parton(-1, x, Q),
		_impl->pdf.parton(-3, x, Q),
	};
}

// Transverse momentum distributions.
FlavorVec ProkudinSfSet::xf1TperpM1(Real x, Real Q_sq) const {
	// Equation [2.A.4].
	return -std::sqrt(E/2.)/(M*SIVERS_M_1)
		*sq(SIVERS_MEAN_K_PERP_SQ)/F1_MEAN_K_PERP_SQ
		*SIVERS_N
		*pow_vec(x, SIVERS_ALPHA)*pow_vec(1. - x, SIVERS_BETA)
		*pow_vec(SIVERS_ALPHA + SIVERS_BETA, SIVERS_ALPHA + SIVERS_BETA)
		*pow_vec(SIVERS_ALPHA, -SIVERS_ALPHA)
		*pow_vec(SIVERS_BETA, -SIVERS_BETA)
		*xf1(x, Q_sq);
}
FlavorVec ProkudinSfSet::xg1(Real x, Real Q_sq) const {
	return x*eval_interp_arr(_impl->interp_g1, x, Q_sq);
}
FlavorVec ProkudinSfSet::xgT(Real x, Real Q_sq) const {
	return eval_interp_arr(_impl->interp_xgT, x, Q_sq);
}
FlavorVec ProkudinSfSet::xh1(Real x, Real Q_sq) const {
	// Use the Soffer bound to get an upper limit on transversity (Equation
	// [2.A.7]).
	return x*H1_N
		*std::pow(x, H1_ALPHA)*std::pow(1. - x, H1_BETA)
		*std::pow(H1_ALPHA + H1_BETA, H1_ALPHA + H1_BETA)
		*std::pow(H1_ALPHA, -H1_ALPHA)
		*std::pow(H1_BETA, -H1_BETA)
		*eval_interp_arr(_impl->interp_sb, x, Q_sq);
}
FlavorVec ProkudinSfSet::xh1M1(Real x, Real Q_sq) const {
	return H1_MEAN_K_PERP_SQ/(2.*sq(M))*xh1(x, Q_sq);
}
FlavorVec ProkudinSfSet::xh1LperpM1(Real x, Real Q_sq) const {
	// Data only exists for up and down quarks.
	return FlavorVec {
		_impl->interp_xh1LperpM1[0]({ x, Q_sq }),
		_impl->interp_xh1LperpM1[1]({ x, Q_sq }),
		0., 0., 0., 0.,
	};
}
FlavorVec ProkudinSfSet::xh1TperpM2(Real x, Real Q_sq) const {
	// Equation [2.A.24].
	return E/(2.*sq(M)*PRETZ_M_TT_SQ)
		*std::pow(PRETZ_MEAN_K_PERP_SQ, 3)/F1_MEAN_K_PERP_SQ
		*PRETZ_N
		*std::pow(x, PRETZ_ALPHA)*std::pow(1. - x, PRETZ_BETA)
		*std::pow(PRETZ_ALPHA + PRETZ_BETA, PRETZ_ALPHA + PRETZ_BETA)
		*std::pow(PRETZ_ALPHA, -PRETZ_ALPHA)
		*std::pow(PRETZ_BETA, -PRETZ_BETA)
		*(xf1(x, Q_sq) - xg1(x, Q_sq));
}
FlavorVec ProkudinSfSet::xh1perpM1(Real x, Real Q_sq) const {
	// Equation [2.A.18].
	return -std::sqrt(E/2.)/(M*BM_M_1)
		*sq(BM_MEAN_K_PERP_SQ)/F1_MEAN_K_PERP_SQ
		*BM_LAMBDA * BM_A
		*pow_vec(x, BM_ALPHA)*pow_vec(1. - x, BM_BETA)
		*pow_vec(BM_ALPHA + BM_BETA, BM_ALPHA + BM_BETA)
		*pow_vec(BM_ALPHA, -BM_ALPHA)
		*pow_vec(BM_BETA, -BM_BETA)
		*xf1(x, Q_sq);
}

