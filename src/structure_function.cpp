#include "sidis/structure_function.hpp"

#include <cmath>

#include <cubature.hpp>

#include "sidis/constant.hpp"
#include "sidis/extra/math.hpp"

// These macros are used to make it a little easier to read the structure
// function definitions below. They compute the convolution integrals of the
// form `C[omega f D]`.
#define CONVOLVE_NUMERIC(weight_type, tmd, ff) \
	convolve_numeric( \
		tmd_set, \
		Weight::weight_type, \
		&TmdSet::x ## tmd, &TmdSet::ff, \
		target, h, x, z, Q_sq, ph_t_sq)
#define CONVOLVE_TILDE_NUMERIC(weight_type, tmd, ff, tmd_tilde, ff_tilde, sign) \
	( \
		(2.*mass(target)*x)/std::sqrt(Q_sq) \
			*CONVOLVE_NUMERIC(weight_type, tmd, ff) \
		+ (sign)*(2.*mass(h))/(z*std::sqrt(Q_sq)) \
			*CONVOLVE_NUMERIC(weight_type, tmd_tilde, ff_tilde))

#define CONVOLVE_GAUSSIAN(weight_type, tmd, ff) \
	convolve_gaussian( \
		tmd_set, \
		Weight::weight_type, \
		tmd_set.mean_ ## tmd.data(), &GaussianTmdSet::x ## tmd, \
		tmd_set.mean_ ## ff.data(), &GaussianTmdSet::ff, \
		target, h, x, z, Q_sq, ph_t_sq)
#define CONVOLVE_TILDE_GAUSSIAN(weight_type, tmd, ff, tmd_tilde, ff_tilde, sign) \
	( \
		(2.*mass(target)*x)/std::sqrt(Q_sq) \
			*CONVOLVE_GAUSSIAN(weight_type, tmd, ff) \
		+ (sign)*(2.*mass(h))/(z*std::sqrt(Q_sq)) \
			*CONVOLVE_GAUSSIAN(weight_type, tmd_tilde, ff_tilde))

#define CONVOLVE(METHOD, weight_type, tmd, ff) \
	(CONVOLVE_##METHOD(weight_type, tmd, ff))
#define CONVOLVE_TILDE(METHOD, weight_type, tmd, ff, tmd_tilde, ff_tilde, sign) \
	(CONVOLVE_TILDE_##METHOD(weight_type, tmd, ff, tmd_tilde, ff_tilde, sign))

// Definitions of structure functions in terms of convolutions on TMDs and FFs.
#define M_F_UUT(METHOD) (CONVOLVE(METHOD, W0, f1, D1))
#define M_F_UU_COS_PHIH(METHOD) ( \
	CONVOLVE_TILDE(METHOD, WA1, h, H1perp, f1, Dperp_tilde, +1) \
	- CONVOLVE_TILDE(METHOD, WB1, fperp, D1, h1perp, H_tilde, +1))
#define M_F_UU_COS_2PHIH(METHOD) (CONVOLVE(METHOD, WAB2, h1perp, H1perp))
#define M_F_UL_SIN_PHIH(METHOD) ( \
	CONVOLVE_TILDE(METHOD, WA1, hL, H1perp, g1, Gperp_tilde, +1) \
	+ CONVOLVE_TILDE(METHOD, WB1, fLperp, D1, h1Lperp, H_tilde, -1))
#define M_F_UL_SIN_2PHIH(METHOD) (CONVOLVE(METHOD, WAB2, h1Lperp, H1perp))
#define M_F_UTT_SIN_PHIH_M_PHIS(METHOD) (-CONVOLVE(METHOD, WB1, f1Tperp, D1))
// This follows equation [2.7.8a], but we choose not to simplify using WW- type
// approximations because the simplification would involve TMDs "before"
// integration, where we have chosen to apply WW-type approximations only
// "after" integration. This means that the following convolution may be less
// efficient than it could be otherwise.
#define M_F_UT_SIN_2PHIH_M_PHIS(METHOD) ( \
	0.5*CONVOLVE_TILDE(METHOD, WAB2, hT, H1perp, g1Tperp, Gperp_tilde, +1) \
	+ 0.5*CONVOLVE_TILDE(METHOD, WAB2, hTperp, H1perp, f1Tperp, Dperp_tilde, -1) \
	+ CONVOLVE_TILDE(METHOD, WC2, fTperp, D1, h1Tperp, H_tilde, -1))
#define M_F_UT_SIN_3PHIH_M_PHIS(METHOD) (CONVOLVE(METHOD, W3, h1Tperp, H1perp))
#define M_F_UT_SIN_PHIS(METHOD) ( \
	CONVOLVE_TILDE(METHOD, W0, fT, D1, h1, H_tilde, -1) \
	- 0.5*CONVOLVE_TILDE(METHOD, WB2, hT, H1perp, g1Tperp, Gperp_tilde, +1) \
	+ 0.5*CONVOLVE_TILDE(METHOD, WB2, hTperp, H1perp, f1Tperp, Dperp_tilde, -1))
#define M_F_UT_SIN_PHIH_P_PHIS(METHOD) (CONVOLVE(METHOD, WA1, h1, H1perp))
#define M_F_LU_SIN_PHIH(METHOD) ( \
	CONVOLVE_TILDE(METHOD, WA1, e, H1perp, f1, Gperp_tilde, +1) \
	+ CONVOLVE_TILDE(METHOD, WB1, gperp, D1, h1perp, E_tilde, +1))
#define M_F_LL(METHOD) (CONVOLVE(METHOD, W0, g1, D1))
#define M_F_LL_COS_PHIH(METHOD) ( \
	-CONVOLVE_TILDE(METHOD, WA1, eL, H1perp, g1, Dperp_tilde, -1) \
	- CONVOLVE_TILDE(METHOD, WB1, gLperp, D1, h1Lperp, E_tilde, +1))
#define M_F_LT_COS_PHIH_M_PHIS(METHOD) (CONVOLVE(METHOD, WB1, g1Tperp, D1))
#define M_F_LT_COS_2PHIH_M_PHIS(METHOD) ( \
	-0.5*CONVOLVE_TILDE(METHOD, WAB2, eT, H1perp, g1Tperp, Dperp_tilde, -1) \
	+ 0.5*CONVOLVE_TILDE(METHOD, WAB2, eTperp, H1perp, f1Tperp, Gperp_tilde, +1) \
	- CONVOLVE_TILDE(METHOD, WC2, gTperp, D1, h1Tperp, E_tilde, +1))
#define M_F_LT_COS_PHIS(METHOD) ( \
	-CONVOLVE_TILDE(METHOD, W0, gT, D1, h1, E_tilde, +1) \
	+ 0.5*CONVOLVE_TILDE(METHOD, WB2, eT, H1perp, g1Tperp, Dperp_tilde, -1) \
	+ 0.5*CONVOLVE_TILDE(METHOD, WB2, eTperp, H1perp, f1Tperp, Gperp_tilde, +1))

// WW approximations for the structure functions.
#define M_F_UUT_WW(METHOD) (CONVOLVE(METHOD, W0, f1, D1))
#define M_F_UU_COS_PHIH_WW(METHOD) ((2.*mass(target)*x)/std::sqrt(Q_sq)*( \
	CONVOLVE(METHOD, WA1, h, H1perp) - CONVOLVE(METHOD, WB1, fperp, D1)))
#define M_F_UU_COS_2PHIH_WW(METHOD) (CONVOLVE(METHOD, WAB2, h1perp, H1perp))
#define M_F_UL_SIN_PHIH_WW(METHOD) ((2.*mass(target)*x)/std::sqrt(Q_sq)* \
	CONVOLVE(METHOD, WA1, hL, H1perp))
#define M_F_UL_SIN_2PHIH_WW(METHOD) (CONVOLVE(METHOD, WAB2, h1Lperp, H1perp))
#define M_F_UTT_SIN_PHIH_M_PHIS_WW(METHOD) (-CONVOLVE(METHOD, WB1, f1Tperp, D1))
#define M_F_UT_SIN_2PHIH_M_PHIS_WW(METHOD) ((2.*mass(target)*x)/std::sqrt(Q_sq)*( \
	0.5*CONVOLVE(METHOD, WAB2, hT, H1perp) \
	+ 0.5*CONVOLVE(METHOD, WAB2, hTperp, H1perp) \
	+ CONVOLVE(METHOD, WC2, fTperp, D1)))
#define M_F_UT_SIN_3PHIH_M_PHIS_WW(METHOD) (CONVOLVE(METHOD, W3, h1Tperp, H1perp))
#define M_F_UT_SIN_PHIS_WW(METHOD) ((2.*mass(target)*x)/std::sqrt(Q_sq)*( \
	CONVOLVE(METHOD, W0, fT, D1) \
	- 0.5*CONVOLVE(METHOD, WB2, hT, H1perp) \
	+ 0.5*CONVOLVE(METHOD, WB2, hTperp, H1perp)))
#define M_F_UT_SIN_PHIH_P_PHIS_WW(METHOD) (CONVOLVE(METHOD, WA1, h1, H1perp))
#define M_F_LU_SIN_PHIH_WW(METHOD) (0.)
#define M_F_LL_WW(METHOD) (CONVOLVE(METHOD, W0, g1, D1))
#define M_F_LL_COS_PHIH_WW(METHOD) (-(2.*mass(target)*x)/std::sqrt(Q_sq)* \
	CONVOLVE(METHOD, WB1, gLperp, D1))
#define M_F_LT_COS_PHIH_M_PHIS_WW(METHOD) (CONVOLVE(METHOD, WB1, g1Tperp, D1))
#define M_F_LT_COS_2PHIH_M_PHIS_WW(METHOD) (-(2.*mass(target)*x)/std::sqrt(Q_sq)* \
	CONVOLVE(METHOD, WC2, gTperp, D1))
#define M_F_LT_COS_PHIS_WW(METHOD) (-(2.*mass(target)*x)/std::sqrt(Q_sq)* \
	CONVOLVE(METHOD, W0, gT, D1))

using namespace sidis;
using namespace sidis::math;
using namespace sidis::sf;

namespace {

using Tmd = Real (TmdSet::*)(unsigned, Real, Real, Real) const;
using Ff = Real (TmdSet::*)(part::Hadron, unsigned, Real, Real, Real) const;
using GaussianTmd = Real (GaussianTmdSet::*)(unsigned, Real, Real) const;
using GaussianFf = Real (GaussianTmdSet::*)(part::Hadron, unsigned, Real, Real) const;

enum class Weight {
	W0,
	WA1,
	WB1,
	WA2,
	WB2,
	WAB2,
	WC2,
	W3,
};

Real convolve_numeric(
		TmdSet const& tmd_set,
		Weight weight_type,
		Tmd tmd, Ff ff,
		part::Nucleus target, part::Hadron h,
		Real x, Real z, Real Q_sq, Real ph_t_sq) {
	Real ph_t = std::sqrt(ph_t_sq);
	Real M = mass(target);
	Real mh = mass(h);
	cubature::EstErr<Real> result = cubature::cubature(
		[&](cubature::Point<2, Real> k_perp_polar) {
			Real dot_k_perp = k_perp_polar[1]*std::cos(k_perp_polar[0]);
			Real dot_p_perp = ph_t - z*dot_k_perp;
			Real k_perp = k_perp_polar[1];
			Real k_perp_sq = sq(k_perp);
			Real p_perp_sq = ph_t_sq + sq(z)*k_perp_sq
				- 2.*ph_t*z*dot_k_perp;
			Real dot_p_k_perp = ph_t*dot_k_perp - z*k_perp_sq;
			Real jacobian = k_perp;
			Real weight;
			switch (weight_type) {
			case Weight::W0:
				weight = 1.;
				break;
			case Weight::WA1:
				weight = dot_p_perp/(z*mh);
				break;
			case Weight::WB1:
				weight = dot_k_perp/M;
				break;
			case Weight::WA2:
				weight = (2.*dot_p_perp*dot_k_perp)/(z*M*mh);
				break;
			case Weight::WB2:
				weight = -dot_p_k_perp/(z*M*mh);
				break;
			case Weight::WAB2:
				weight = (2.*dot_p_perp*dot_k_perp - dot_p_k_perp)/(z*M*mh);
				break;
			case Weight::WC2:
				weight = (2.*sq(dot_k_perp) - k_perp_sq)/(2.*sq(M));
				break;
			case Weight::W3:
				weight = (
					4.*dot_p_perp*sq(dot_k_perp)
					- dot_k_perp*dot_p_k_perp
					- dot_p_perp*k_perp_sq)/(2.*z*sq(M)*mh);
				break;
			default:
				// Unknown weight.
				weight = 0.;
			}
			Real integrand = 0.;
			for (unsigned fl = 0; fl < tmd_set.flavor_count; ++fl) {
				integrand += sq(tmd_set.charge(fl))
					*(tmd_set.*tmd)(fl, x, Q_sq, sq(k_perp))
					*(tmd_set.*ff)(h, fl, z, Q_sq, p_perp_sq);
			}
			return jacobian*weight*integrand;
		},
		cubature::Point<2, Real>{ 0., 0. },
		// TODO: Choose the upper bound on `k_perp` using a better method.
		cubature::Point<2, Real>{ 2. * PI, 10. },
		100000, 0., 1e-6);
	return result.val;
}

Real convolve_gaussian(
		GaussianTmdSet const& tmd_set,
		Weight weight_type,
		Real const* mean_tmd, GaussianTmd tmd,
		Real const* mean_ff, GaussianFf ff,
		part::Nucleus target, part::Hadron h,
		Real x, Real z, Real Q_sq, Real ph_t_sq) {
	Real M = mass(target);
	Real mh = mass(h);
	Real ph_t = std::sqrt(ph_t_sq);
	Real result = 0.;
	for (unsigned fl = 0; fl < tmd_set.flavor_count; ++fl) {
		Real mean = mean_ff[fl] + sq(z)*mean_tmd[fl];
		if (std::isinf(mean)) {
			continue;
		}
		Real gaussian = std::exp(-ph_t_sq/mean)/(PI*mean);
		// Analytically evaluate the convolution integral for the TMD and FF
		// Gaussians with the given means.
		Real weight;
		switch (weight_type) {
		case Weight::W0:
			weight = 1.;
			break;
		case Weight::WA1:
			weight = (mean_ff[fl]*ph_t)/(mh*mean*z);
			break;
		case Weight::WB1:
			weight = (mean_tmd[fl]*ph_t*z)/(M*mean);
			break;
		case Weight::WA2:
			weight = (mean_tmd[fl]*mean_ff[fl]*(2.*ph_t_sq - mean))/(M*mh*sq(mean));
			break;
		case Weight::WB2:
			weight = (mean_tmd[fl]*mean_ff[fl]*(mean - ph_t_sq))/(M*mh*sq(mean));
			break;
		case Weight::WAB2:
			weight = (mean_tmd[fl]*mean_ff[fl]*ph_t_sq)/(M*mh*sq(mean));
			break;
		case Weight::WC2:
			weight = (sq(mean_tmd[fl])*ph_t_sq*sq(z))/(2.*sq(M)*sq(mean));
			break;
		case Weight::W3:
			weight = (sq(mean_tmd[fl])*mean_ff[fl]*ph_t*ph_t_sq*z)
				/(2.*sq(M)*mh*mean*sq(mean));
			break;
		default:
			// Unknown integrand.
			weight = 0.;
		}
		result += sq(tmd_set.charge(fl))
			*weight*gaussian
			*(tmd_set.*tmd)(fl, x, Q_sq)
			*(tmd_set.*ff)(h, fl, z, Q_sq);
	}
	return result;
}

}

Real SfSet::F_UUL(part::Hadron, Real, Real, Real, Real) const {
	return 0.;
}
Real SfSet::F_UUT(part::Hadron, Real, Real, Real, Real) const {
	return 0.;
}
Real SfSet::F_UU_cos_phih(part::Hadron, Real, Real, Real, Real) const {
	return 0.;
}
Real SfSet::F_UU_cos_2phih(part::Hadron, Real, Real, Real, Real) const {
	return 0.;
}

Real SfSet::F_UL_sin_phih(part::Hadron, Real, Real, Real, Real) const {
	return 0.;
}
Real SfSet::F_UL_sin_2phih(part::Hadron, Real, Real, Real, Real) const {
	return 0.;
}

Real SfSet::F_UTL_sin_phih_m_phis(part::Hadron, Real, Real, Real, Real) const {
	return 0.;
}
Real SfSet::F_UTT_sin_phih_m_phis(part::Hadron, Real, Real, Real, Real) const {
	return 0.;
}
Real SfSet::F_UT_sin_2phih_m_phis(part::Hadron, Real, Real, Real, Real) const {
	return 0.;
}
Real SfSet::F_UT_sin_3phih_m_phis(part::Hadron, Real, Real, Real, Real) const {
	return 0.;
}
Real SfSet::F_UT_sin_phis(part::Hadron, Real, Real, Real, Real) const {
	return 0.;
}
Real SfSet::F_UT_sin_phih_p_phis(part::Hadron, Real, Real, Real, Real) const {
	return 0.;
}

Real SfSet::F_LU_sin_phih(part::Hadron, Real, Real, Real, Real) const {
	return 0.;
}

Real SfSet::F_LL(part::Hadron, Real, Real, Real, Real) const {
	return 0.;
}
Real SfSet::F_LL_cos_phih(part::Hadron, Real, Real, Real, Real) const {
	return 0.;
}

Real SfSet::F_LT_cos_phih_m_phis(part::Hadron, Real, Real, Real, Real) const {
	return 0.;
}
Real SfSet::F_LT_cos_2phih_m_phis(part::Hadron, Real, Real, Real, Real) const {
	return 0.;
}
Real SfSet::F_LT_cos_phis(part::Hadron, Real, Real, Real, Real) const {
	return 0.;
}

SfBaseUU SfSet::sf_base_uu(part::Hadron h, Real x, Real z, Real Q_sq, Real ph_t_sq) const {
	return {
		F_UUL(h, x, z, Q_sq, ph_t_sq),
		F_UUT(h, x, z, Q_sq, ph_t_sq),
		F_UU_cos_phih(h, x, z, Q_sq, ph_t_sq),
		F_UU_cos_2phih(h, x, z, Q_sq, ph_t_sq),
	};
}
SfBaseUL SfSet::sf_base_ul(part::Hadron h, Real x, Real z, Real Q_sq, Real ph_t_sq) const {
	return {
		F_UL_sin_phih(h, x, z, Q_sq, ph_t_sq),
		F_UL_sin_2phih(h, x, z, Q_sq, ph_t_sq),
	};
}
SfBaseUT SfSet::sf_base_ut(part::Hadron h, Real x, Real z, Real Q_sq, Real ph_t_sq) const {
	return {
		F_UTL_sin_phih_m_phis(h, x, z, Q_sq, ph_t_sq),
		F_UTT_sin_phih_m_phis(h, x, z, Q_sq, ph_t_sq),
		F_UT_sin_2phih_m_phis(h, x, z, Q_sq, ph_t_sq),
		F_UT_sin_3phih_m_phis(h, x, z, Q_sq, ph_t_sq),
		F_UT_sin_phis(h, x, z, Q_sq, ph_t_sq),
		F_UT_sin_phih_p_phis(h, x, z, Q_sq, ph_t_sq),
	};
}
SfBaseUP SfSet::sf_base_up(part::Hadron h, Real x, Real z, Real Q_sq, Real ph_t_sq) const {
	return {
		sf_base_ul(h, x, z, Q_sq, ph_t_sq),
		sf_base_ut(h, x, z, Q_sq, ph_t_sq),
	};
}
SfBaseLU SfSet::sf_base_lu(part::Hadron h, Real x, Real z, Real Q_sq, Real ph_t_sq) const {
	return {
		F_LU_sin_phih(h, x, z, Q_sq, ph_t_sq),
	};
}
SfBaseLL SfSet::sf_base_ll(part::Hadron h, Real x, Real z, Real Q_sq, Real ph_t_sq) const {
	return {
		F_LL(h, x, z, Q_sq, ph_t_sq),
		F_LL_cos_phih(h, x, z, Q_sq, ph_t_sq),
	};
}
SfBaseLT SfSet::sf_base_lt(part::Hadron h, Real x, Real z, Real Q_sq, Real ph_t_sq) const {
	return {
		F_LT_cos_phih_m_phis(h, x, z, Q_sq, ph_t_sq),
		F_LT_cos_2phih_m_phis(h, x, z, Q_sq, ph_t_sq),
		F_LT_cos_phis(h, x, z, Q_sq, ph_t_sq),
	};
}
SfBaseLP SfSet::sf_base_lp(part::Hadron h, Real x, Real z, Real Q_sq, Real ph_t_sq) const {
	return {
		sf_base_ll(h, x, z, Q_sq, ph_t_sq),
		sf_base_lt(h, x, z, Q_sq, ph_t_sq),
	};
}

SfUU SfSet::sf_uu(part::Hadron h, Real x, Real z, Real Q_sq, Real ph_t_sq) const {
	return {
		sf_base_uu(h, x, z, Q_sq, ph_t_sq),
	};
}
SfUL SfSet::sf_ul(part::Hadron h, Real x, Real z, Real Q_sq, Real ph_t_sq) const {
	return {
		sf_base_uu(h, x, z, Q_sq, ph_t_sq),
		sf_base_ul(h, x, z, Q_sq, ph_t_sq),
	};
}
SfUT SfSet::sf_ut(part::Hadron h, Real x, Real z, Real Q_sq, Real ph_t_sq) const {
	return {
		sf_base_uu(h, x, z, Q_sq, ph_t_sq),
		sf_base_ut(h, x, z, Q_sq, ph_t_sq),
	};
}
SfUP SfSet::sf_up(part::Hadron h, Real x, Real z, Real Q_sq, Real ph_t_sq) const {
	return {
		sf_base_uu(h, x, z, Q_sq, ph_t_sq),
		sf_base_ul(h, x, z, Q_sq, ph_t_sq),
		sf_base_ut(h, x, z, Q_sq, ph_t_sq),
	};
}
SfLU SfSet::sf_lu(part::Hadron h, Real x, Real z, Real Q_sq, Real ph_t_sq) const {
	return {
		sf_base_uu(h, x, z, Q_sq, ph_t_sq),
		sf_base_lu(h, x, z, Q_sq, ph_t_sq),
	};
}
SfLL SfSet::sf_ll(part::Hadron h, Real x, Real z, Real Q_sq, Real ph_t_sq) const {
	return {
		sf_base_uu(h, x, z, Q_sq, ph_t_sq),
		sf_base_ul(h, x, z, Q_sq, ph_t_sq),
		sf_base_lu(h, x, z, Q_sq, ph_t_sq),
		sf_base_ll(h, x, z, Q_sq, ph_t_sq),
	};
}
SfLT SfSet::sf_lt(part::Hadron h, Real x, Real z, Real Q_sq, Real ph_t_sq) const {
	return {
		sf_base_uu(h, x, z, Q_sq, ph_t_sq),
		sf_base_ut(h, x, z, Q_sq, ph_t_sq),
		sf_base_lu(h, x, z, Q_sq, ph_t_sq),
		sf_base_lt(h, x, z, Q_sq, ph_t_sq),
	};
}
SfLP SfSet::sf_lp(part::Hadron h, Real x, Real z, Real Q_sq, Real ph_t_sq) const {
	return {
		sf_base_uu(h, x, z, Q_sq, ph_t_sq),
		sf_base_ul(h, x, z, Q_sq, ph_t_sq),
		sf_base_ut(h, x, z, Q_sq, ph_t_sq),
		sf_base_lu(h, x, z, Q_sq, ph_t_sq),
		sf_base_ll(h, x, z, Q_sq, ph_t_sq),
		sf_base_lt(h, x, z, Q_sq, ph_t_sq),
	};
}

// Full structure function calculations from equations [2.17], [2.18].
// TODO: Find expression for longitudinally-polarized photon terms.
Real TmdSfSet::F_UUL(part::Hadron, Real, Real, Real, Real) const {
	return 0.;
}
Real TmdSfSet::F_UUT(part::Hadron h, Real x, Real z, Real Q_sq, Real ph_t_sq) const {
	return M_F_UUT(NUMERIC);
}
Real TmdSfSet::F_UU_cos_phih(part::Hadron h, Real x, Real z, Real Q_sq, Real ph_t_sq) const {
	return M_F_UU_COS_PHIH(NUMERIC);
}
Real TmdSfSet::F_UU_cos_2phih(part::Hadron h, Real x, Real z, Real Q_sq, Real ph_t_sq) const {
	return M_F_UU_COS_2PHIH(NUMERIC);
}

Real TmdSfSet::F_UL_sin_phih(part::Hadron h, Real x, Real z, Real Q_sq, Real ph_t_sq) const {
	return M_F_UL_SIN_PHIH(NUMERIC);
}
Real TmdSfSet::F_UL_sin_2phih(part::Hadron h, Real x, Real z, Real Q_sq, Real ph_t_sq) const {
	return M_F_UL_SIN_2PHIH(NUMERIC);
}

Real TmdSfSet::F_UTL_sin_phih_m_phis(part::Hadron, Real, Real, Real, Real) const {
	return 0.;
}
Real TmdSfSet::F_UTT_sin_phih_m_phis(part::Hadron h, Real x, Real z, Real Q_sq, Real ph_t_sq) const {
	return M_F_UTT_SIN_PHIH_M_PHIS(NUMERIC);
}
Real TmdSfSet::F_UT_sin_2phih_m_phis(part::Hadron h, Real x, Real z, Real Q_sq, Real ph_t_sq) const {
	return M_F_UT_SIN_2PHIH_M_PHIS(NUMERIC);
}
Real TmdSfSet::F_UT_sin_3phih_m_phis(part::Hadron h, Real x, Real z, Real Q_sq, Real ph_t_sq) const {
	return M_F_UT_SIN_3PHIH_M_PHIS(NUMERIC);
}
Real TmdSfSet::F_UT_sin_phis(part::Hadron h, Real x, Real z, Real Q_sq, Real ph_t_sq) const {
	return M_F_UT_SIN_PHIS(NUMERIC);
}
Real TmdSfSet::F_UT_sin_phih_p_phis(part::Hadron h, Real x, Real z, Real Q_sq, Real ph_t_sq) const {
	return M_F_UT_SIN_PHIH_P_PHIS(NUMERIC);
}

Real TmdSfSet::F_LU_sin_phih(part::Hadron h, Real x, Real z, Real Q_sq, Real ph_t_sq) const {
	return M_F_LU_SIN_PHIH(NUMERIC);
}

Real TmdSfSet::F_LL(part::Hadron h, Real x, Real z, Real Q_sq, Real ph_t_sq) const {
	return M_F_LL(NUMERIC);
}
Real TmdSfSet::F_LL_cos_phih(part::Hadron h, Real x, Real z, Real Q_sq, Real ph_t_sq) const {
	return M_F_LL_COS_PHIH(NUMERIC);
}

Real TmdSfSet::F_LT_cos_phih_m_phis(part::Hadron h, Real x, Real z, Real Q_sq, Real ph_t_sq) const {
	return M_F_LT_COS_PHIH_M_PHIS(NUMERIC);
}
Real TmdSfSet::F_LT_cos_2phih_m_phis(part::Hadron h, Real x, Real z, Real Q_sq, Real ph_t_sq) const {
	return M_F_LT_COS_2PHIH_M_PHIS(NUMERIC);
}
Real TmdSfSet::F_LT_cos_phis(part::Hadron h, Real x, Real z, Real Q_sq, Real ph_t_sq) const {
	return M_F_LT_COS_PHIS(NUMERIC);
}

// Gaussian approximation.
Real GaussianTmdSfSet::F_UUL(part::Hadron, Real, Real, Real, Real) const {
	return 0.;
}
Real GaussianTmdSfSet::F_UUT(part::Hadron h, Real x, Real z, Real Q_sq, Real ph_t_sq) const {
	return M_F_UUT(GAUSSIAN);
}
Real GaussianTmdSfSet::F_UU_cos_phih(part::Hadron h, Real x, Real z, Real Q_sq, Real ph_t_sq) const {
	return M_F_UU_COS_PHIH(GAUSSIAN);
}
Real GaussianTmdSfSet::F_UU_cos_2phih(part::Hadron h, Real x, Real z, Real Q_sq, Real ph_t_sq) const {
	return M_F_UU_COS_2PHIH(GAUSSIAN);
}

Real GaussianTmdSfSet::F_UL_sin_phih(part::Hadron h, Real x, Real z, Real Q_sq, Real ph_t_sq) const {
	return M_F_UL_SIN_PHIH(GAUSSIAN);
}
Real GaussianTmdSfSet::F_UL_sin_2phih(part::Hadron h, Real x, Real z, Real Q_sq, Real ph_t_sq) const {
	return M_F_UL_SIN_2PHIH(GAUSSIAN);
}

Real GaussianTmdSfSet::F_UTL_sin_phih_m_phis(part::Hadron, Real, Real, Real, Real) const {
	return 0.;
}
Real GaussianTmdSfSet::F_UTT_sin_phih_m_phis(part::Hadron h, Real x, Real z, Real Q_sq, Real ph_t_sq) const {
	return M_F_UTT_SIN_PHIH_M_PHIS(GAUSSIAN);
}
Real GaussianTmdSfSet::F_UT_sin_2phih_m_phis(part::Hadron h, Real x, Real z, Real Q_sq, Real ph_t_sq) const {
	return M_F_UT_SIN_2PHIH_M_PHIS(GAUSSIAN);
}
Real GaussianTmdSfSet::F_UT_sin_3phih_m_phis(part::Hadron h, Real x, Real z, Real Q_sq, Real ph_t_sq) const {
	return M_F_UT_SIN_3PHIH_M_PHIS(GAUSSIAN);
}
Real GaussianTmdSfSet::F_UT_sin_phis(part::Hadron h, Real x, Real z, Real Q_sq, Real ph_t_sq) const {
	return M_F_UT_SIN_PHIS(GAUSSIAN);
}
Real GaussianTmdSfSet::F_UT_sin_phih_p_phis(part::Hadron h, Real x, Real z, Real Q_sq, Real ph_t_sq) const {
	return M_F_UT_SIN_PHIH_P_PHIS(GAUSSIAN);
}

Real GaussianTmdSfSet::F_LU_sin_phih(part::Hadron h, Real x, Real z, Real Q_sq, Real ph_t_sq) const {
	return M_F_LU_SIN_PHIH(GAUSSIAN);
}

Real GaussianTmdSfSet::F_LL(part::Hadron h, Real x, Real z, Real Q_sq, Real ph_t_sq) const {
	return M_F_LL(GAUSSIAN);
}
Real GaussianTmdSfSet::F_LL_cos_phih(part::Hadron h, Real x, Real z, Real Q_sq, Real ph_t_sq) const {
	return M_F_LL_COS_PHIH(GAUSSIAN);
}

Real GaussianTmdSfSet::F_LT_cos_phih_m_phis(part::Hadron h, Real x, Real z, Real Q_sq, Real ph_t_sq) const {
	return M_F_LT_COS_PHIH_M_PHIS(GAUSSIAN);
}
Real GaussianTmdSfSet::F_LT_cos_2phih_m_phis(part::Hadron h, Real x, Real z, Real Q_sq, Real ph_t_sq) const {
	return M_F_LT_COS_2PHIH_M_PHIS(GAUSSIAN);
}
Real GaussianTmdSfSet::F_LT_cos_phis(part::Hadron h, Real x, Real z, Real Q_sq, Real ph_t_sq) const {
	return M_F_LT_COS_PHIS(GAUSSIAN);
}

// WW-type approximation.
Real WwTmdSfSet::F_UUL(part::Hadron, Real, Real, Real, Real) const {
	return 0.;
}
Real WwTmdSfSet::F_UUT(part::Hadron h, Real x, Real z, Real Q_sq, Real ph_t_sq) const {
	return M_F_UUT_WW(NUMERIC);
}
Real WwTmdSfSet::F_UU_cos_phih(part::Hadron h, Real x, Real z, Real Q_sq, Real ph_t_sq) const {
	return M_F_UU_COS_PHIH_WW(NUMERIC);
}
Real WwTmdSfSet::F_UU_cos_2phih(part::Hadron h, Real x, Real z, Real Q_sq, Real ph_t_sq) const {
	return M_F_UU_COS_2PHIH_WW(NUMERIC);
}

Real WwTmdSfSet::F_UL_sin_phih(part::Hadron h, Real x, Real z, Real Q_sq, Real ph_t_sq) const {
	return M_F_UL_SIN_PHIH_WW(NUMERIC);
}
Real WwTmdSfSet::F_UL_sin_2phih(part::Hadron h, Real x, Real z, Real Q_sq, Real ph_t_sq) const {
	return M_F_UL_SIN_2PHIH_WW(NUMERIC);
}

Real WwTmdSfSet::F_UTL_sin_phih_m_phis(part::Hadron, Real, Real, Real, Real) const {
	return 0.;
}
Real WwTmdSfSet::F_UTT_sin_phih_m_phis(part::Hadron h, Real x, Real z, Real Q_sq, Real ph_t_sq) const {
	return M_F_UTT_SIN_PHIH_M_PHIS_WW(NUMERIC);
}
Real WwTmdSfSet::F_UT_sin_2phih_m_phis(part::Hadron h, Real x, Real z, Real Q_sq, Real ph_t_sq) const {
	return M_F_UT_SIN_2PHIH_M_PHIS_WW(NUMERIC);
}
Real WwTmdSfSet::F_UT_sin_3phih_m_phis(part::Hadron h, Real x, Real z, Real Q_sq, Real ph_t_sq) const {
	return M_F_UT_SIN_3PHIH_M_PHIS_WW(NUMERIC);
}
Real WwTmdSfSet::F_UT_sin_phis(part::Hadron h, Real x, Real z, Real Q_sq, Real ph_t_sq) const {
	return M_F_UT_SIN_PHIS_WW(NUMERIC);
}
Real WwTmdSfSet::F_UT_sin_phih_p_phis(part::Hadron h, Real x, Real z, Real Q_sq, Real ph_t_sq) const {
	return M_F_UT_SIN_PHIH_P_PHIS_WW(NUMERIC);
}

Real WwTmdSfSet::F_LU_sin_phih(part::Hadron, Real, Real, Real, Real) const {
	return M_F_LU_SIN_PHIH_WW(NUMERIC);
}

Real WwTmdSfSet::F_LL(part::Hadron h, Real x, Real z, Real Q_sq, Real ph_t_sq) const {
	return M_F_LL_WW(NUMERIC);
}
Real WwTmdSfSet::F_LL_cos_phih(part::Hadron h, Real x, Real z, Real Q_sq, Real ph_t_sq) const {
	return M_F_LL_COS_PHIH_WW(NUMERIC);
}

Real WwTmdSfSet::F_LT_cos_phih_m_phis(part::Hadron h, Real x, Real z, Real Q_sq, Real ph_t_sq) const {
	return M_F_LT_COS_PHIH_M_PHIS_WW(NUMERIC);
}
Real WwTmdSfSet::F_LT_cos_2phih_m_phis(part::Hadron h, Real x, Real z, Real Q_sq, Real ph_t_sq) const {
	return M_F_LT_COS_2PHIH_M_PHIS_WW(NUMERIC);
}
Real WwTmdSfSet::F_LT_cos_phis(part::Hadron h, Real x, Real z, Real Q_sq, Real ph_t_sq) const {
	return M_F_LT_COS_PHIS_WW(NUMERIC);
}

// WW-type approximation combined with Gaussian TMDs.
Real GaussianWwTmdSfSet::F_UUL(part::Hadron, Real, Real, Real, Real) const {
	return 0.;
}
Real GaussianWwTmdSfSet::F_UUT(part::Hadron h, Real x, Real z, Real Q_sq, Real ph_t_sq) const {
	return M_F_UUT_WW(GAUSSIAN);
}
Real GaussianWwTmdSfSet::F_UU_cos_phih(part::Hadron h, Real x, Real z, Real Q_sq, Real ph_t_sq) const {
	return M_F_UU_COS_PHIH_WW(GAUSSIAN);
}
Real GaussianWwTmdSfSet::F_UU_cos_2phih(part::Hadron h, Real x, Real z, Real Q_sq, Real ph_t_sq) const {
	return M_F_UU_COS_2PHIH_WW(GAUSSIAN);
}

Real GaussianWwTmdSfSet::F_UL_sin_phih(part::Hadron h, Real x, Real z, Real Q_sq, Real ph_t_sq) const {
	return M_F_UL_SIN_PHIH_WW(GAUSSIAN);
}
Real GaussianWwTmdSfSet::F_UL_sin_2phih(part::Hadron h, Real x, Real z, Real Q_sq, Real ph_t_sq) const {
	return M_F_UL_SIN_2PHIH_WW(GAUSSIAN);
}

Real GaussianWwTmdSfSet::F_UTL_sin_phih_m_phis(part::Hadron, Real, Real, Real, Real) const {
	return 0.;
}
Real GaussianWwTmdSfSet::F_UTT_sin_phih_m_phis(part::Hadron h, Real x, Real z, Real Q_sq, Real ph_t_sq) const {
	return M_F_UTT_SIN_PHIH_M_PHIS_WW(GAUSSIAN);
}
Real GaussianWwTmdSfSet::F_UT_sin_2phih_m_phis(part::Hadron h, Real x, Real z, Real Q_sq, Real ph_t_sq) const {
	return M_F_UT_SIN_2PHIH_M_PHIS_WW(GAUSSIAN);
}
Real GaussianWwTmdSfSet::F_UT_sin_3phih_m_phis(part::Hadron h, Real x, Real z, Real Q_sq, Real ph_t_sq) const {
	return M_F_UT_SIN_3PHIH_M_PHIS_WW(GAUSSIAN);
}
Real GaussianWwTmdSfSet::F_UT_sin_phis(part::Hadron h, Real x, Real z, Real Q_sq, Real ph_t_sq) const {
	return M_F_UT_SIN_PHIS_WW(GAUSSIAN);
}
Real GaussianWwTmdSfSet::F_UT_sin_phih_p_phis(part::Hadron h, Real x, Real z, Real Q_sq, Real ph_t_sq) const {
	return M_F_UT_SIN_PHIH_P_PHIS_WW(GAUSSIAN);
}

Real GaussianWwTmdSfSet::F_LU_sin_phih(part::Hadron, Real, Real, Real, Real) const {
	return M_F_LU_SIN_PHIH_WW(GAUSSIAN);
}

Real GaussianWwTmdSfSet::F_LL(part::Hadron h, Real x, Real z, Real Q_sq, Real ph_t_sq) const {
	return M_F_LL_WW(GAUSSIAN);
}
Real GaussianWwTmdSfSet::F_LL_cos_phih(part::Hadron h, Real x, Real z, Real Q_sq, Real ph_t_sq) const {
	return M_F_LL_COS_PHIH_WW(GAUSSIAN);
}

Real GaussianWwTmdSfSet::F_LT_cos_phih_m_phis(part::Hadron h, Real x, Real z, Real Q_sq, Real ph_t_sq) const {
	return M_F_LT_COS_PHIH_M_PHIS_WW(GAUSSIAN);
}
Real GaussianWwTmdSfSet::F_LT_cos_2phih_m_phis(part::Hadron h, Real x, Real z, Real Q_sq, Real ph_t_sq) const {
	return M_F_LT_COS_2PHIH_M_PHIS_WW(GAUSSIAN);
}
Real GaussianWwTmdSfSet::F_LT_cos_phis(part::Hadron h, Real x, Real z, Real Q_sq, Real ph_t_sq) const {
	return M_F_LT_COS_PHIS_WW(GAUSSIAN);
}

