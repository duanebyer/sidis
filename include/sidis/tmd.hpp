#ifndef SIDIS_Tmd_HPP
#define SIDIS_Tmd_HPP

#include "sidis/constant.hpp"

namespace sidis {
namespace sf {

class TmdSet {
public:
	unsigned const flavor_count;
	constant::Nucleus const target;

	TmdSet(unsigned flavor_count, constant::Nucleus target) :
		flavor_count(flavor_count),
		target(target) { }
	virtual ~TmdSet() = default;

	// Charge of a given flavor (used as a weighting when computing structure
	// functions).
	virtual Real charge(unsigned fl) const = 0;

	// Transverse momentum distribtions (TMDs).
	virtual Real xf1(unsigned fl, Real x, Real Q_sq, Real k_perp_sq) const;
	virtual Real xf1Tperp(unsigned fl, Real x, Real Q_sq, Real k_perp_sq) const;
	virtual Real xfT(unsigned fl, Real x, Real Q_sq, Real k_perp_sq) const;
	virtual Real xfperp(unsigned fl, Real x, Real Q_sq, Real k_perp_sq) const;
	virtual Real xfLperp(unsigned fl, Real x, Real Q_sq, Real k_perp_sq) const;
	virtual Real xfTperp(unsigned fl, Real x, Real Q_sq, Real k_perp_sq) const;
	virtual Real xg1(unsigned fl, Real x, Real Q_sq, Real k_perp_sq) const;
	virtual Real xg1Tperp(unsigned fl, Real x, Real Q_sq, Real k_perp_sq) const;
	virtual Real xgT(unsigned fl, Real x, Real Q_sq, Real k_perp_sq) const;
	virtual Real xgperp(unsigned fl, Real x, Real Q_sq, Real k_perp_sq) const;
	virtual Real xgLperp(unsigned fl, Real x, Real Q_sq, Real k_perp_sq) const;
	virtual Real xgTperp(unsigned fl, Real x, Real Q_sq, Real k_perp_sq) const;
	virtual Real xh1(unsigned fl, Real x, Real Q_sq, Real k_perp_sq) const;
	virtual Real xh1perp(unsigned fl, Real x, Real Q_sq, Real k_perp_sq) const;
	virtual Real xh1Lperp(unsigned fl, Real x, Real Q_sq, Real k_perp_sq) const;
	virtual Real xh1Tperp(unsigned fl, Real x, Real Q_sq, Real k_perp_sq) const;
	virtual Real xh(unsigned fl, Real x, Real Q_sq, Real k_perp_sq) const;
	virtual Real xhL(unsigned fl, Real x, Real Q_sq, Real k_perp_sq) const;
	virtual Real xhT(unsigned fl, Real x, Real Q_sq, Real k_perp_sq) const;
	virtual Real xhTperp(unsigned fl, Real x, Real Q_sq, Real k_perp_sq) const;
	virtual Real xe(unsigned fl, Real x, Real Q_sq, Real k_perp_sq) const;
	virtual Real xeL(unsigned fl, Real x, Real Q_sq, Real k_perp_sq) const;
	virtual Real xeT(unsigned fl, Real x, Real Q_sq, Real k_perp_sq) const;
	virtual Real xeTperp(unsigned fl, Real x, Real Q_sq, Real k_perp_sq) const;

	// Fragmentation functions (FFs).
	virtual Real D1(constant::Hadron h, unsigned fl, Real z, Real Q_sq, Real p_perp_sq) const;
	virtual Real H1perp(constant::Hadron h, unsigned fl, Real z, Real Q_sq, Real p_perp_sq) const;
	virtual Real Dperp_tilde(constant::Hadron h, unsigned fl, Real z, Real Q_sq, Real p_perp_sq) const;
	virtual Real H_tilde(constant::Hadron h, unsigned fl, Real z, Real Q_sq, Real p_perp_sq) const;
	virtual Real Gperp_tilde(constant::Hadron h, unsigned fl, Real z, Real Q_sq, Real p_perp_sq) const;
	virtual Real E_tilde(constant::Hadron h, unsigned fl, Real z, Real Q_sq, Real p_perp_sq) const;
};

class TmdGaussianSet : public TmdSet {
public:
	// TODO: For now, the widths of the Gaussian approximations are required to
	// be independent of flavor, but support for flavor-dependent widths may be
	// added later.
	Real const mean_f1;
	Real const mean_f1Tperp;
	Real const mean_fT;
	Real const mean_fperp;
	Real const mean_fLperp;
	Real const mean_fTperp;
	Real const mean_g1;
	Real const mean_g1Tperp;
	Real const mean_gT;
	Real const mean_gperp;
	Real const mean_gLperp;
	Real const mean_gTperp;
	Real const mean_h1;
	Real const mean_h1perp;
	Real const mean_h1Lperp;
	Real const mean_h1Tperp;
	Real const mean_h;
	Real const mean_hL;
	Real const mean_hT;
	Real const mean_hTperp;
	Real const mean_e;
	Real const mean_eL;
	Real const mean_eT;
	Real const mean_eTperp;

	Real const mean_D1;
	Real const mean_H1perp;
	Real const mean_Dperp_tilde;
	Real const mean_H_tilde;
	Real const mean_Gperp_tilde;
	Real const mean_E_tilde;

	TmdGaussianSet(
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
		Real mean_E_tilde);
	virtual ~TmdGaussianSet() = default;

	// Reduced TMDs (without `k_perp` dependence).
	virtual Real xf1(unsigned fl, Real x, Real Q_sq) const;
	virtual Real xf1Tperp(unsigned fl, Real x, Real Q_sq) const;
	virtual Real xfT(unsigned fl, Real x, Real Q_sq) const;
	virtual Real xfperp(unsigned fl, Real x, Real Q_sq) const;
	virtual Real xfLperp(unsigned fl, Real x, Real Q_sq) const;
	virtual Real xfTperp(unsigned fl, Real x, Real Q_sq) const;
	virtual Real xg1(unsigned fl, Real x, Real Q_sq) const;
	virtual Real xg1Tperp(unsigned fl, Real x, Real Q_sq) const;
	virtual Real xgT(unsigned fl, Real x, Real Q_sq) const;
	virtual Real xgperp(unsigned fl, Real x, Real Q_sq) const;
	virtual Real xgLperp(unsigned fl, Real x, Real Q_sq) const;
	virtual Real xgTperp(unsigned fl, Real x, Real Q_sq) const;
	virtual Real xh1(unsigned fl, Real x, Real Q_sq) const;
	virtual Real xh1perp(unsigned fl, Real x, Real Q_sq) const;
	virtual Real xh1Lperp(unsigned fl, Real x, Real Q_sq) const;
	virtual Real xh1Tperp(unsigned fl, Real x, Real Q_sq) const;
	virtual Real xh(unsigned fl, Real x, Real Q_sq) const;
	virtual Real xhL(unsigned fl, Real x, Real Q_sq) const;
	virtual Real xhT(unsigned fl, Real x, Real Q_sq) const;
	virtual Real xhTperp(unsigned fl, Real x, Real Q_sq) const;
	virtual Real xe(unsigned fl, Real x, Real Q_sq) const;
	virtual Real xeL(unsigned fl, Real x, Real Q_sq) const;
	virtual Real xeT(unsigned fl, Real x, Real Q_sq) const;
	virtual Real xeTperp(unsigned fl, Real x, Real Q_sq) const;

	// Reduced FFs (without `p_perp` dependence).
	virtual Real D1(constant::Hadron h, unsigned fl, Real z, Real Q_sq) const;
	virtual Real H1perp(constant::Hadron h, unsigned fl, Real z, Real Q_sq) const;
	virtual Real Dperp_tilde(constant::Hadron h, unsigned fl, Real z, Real Q_sq) const;
	virtual Real H_tilde(constant::Hadron h, unsigned fl, Real z, Real Q_sq) const;
	virtual Real Gperp_tilde(constant::Hadron h, unsigned fl, Real z, Real Q_sq) const;
	virtual Real E_tilde(constant::Hadron h, unsigned fl, Real z, Real Q_sq) const;

	Real xf1(unsigned fl, Real x, Real Q_sq, Real k_perp_sq) const override;
	Real xf1Tperp(unsigned fl, Real x, Real Q_sq, Real k_perp_sq) const override;
	Real xfT(unsigned fl, Real x, Real Q_sq, Real k_perp_sq) const override;
	Real xfperp(unsigned fl, Real x, Real Q_sq, Real k_perp_sq) const override;
	Real xfLperp(unsigned fl, Real x, Real Q_sq, Real k_perp_sq) const override;
	Real xfTperp(unsigned fl, Real x, Real Q_sq, Real k_perp_sq) const override;
	Real xg1(unsigned fl, Real x, Real Q_sq, Real k_perp_sq) const override;
	Real xg1Tperp(unsigned fl, Real x, Real Q_sq, Real k_perp_sq) const override;
	Real xgT(unsigned fl, Real x, Real Q_sq, Real k_perp_sq) const override;
	Real xgperp(unsigned fl, Real x, Real Q_sq, Real k_perp_sq) const override;
	Real xgLperp(unsigned fl, Real x, Real Q_sq, Real k_perp_sq) const override;
	Real xgTperp(unsigned fl, Real x, Real Q_sq, Real k_perp_sq) const override;
	Real xh1(unsigned fl, Real x, Real Q_sq, Real k_perp_sq) const override;
	Real xh1perp(unsigned fl, Real x, Real Q_sq, Real k_perp_sq) const override;
	Real xh1Lperp(unsigned fl, Real x, Real Q_sq, Real k_perp_sq) const override;
	Real xh1Tperp(unsigned fl, Real x, Real Q_sq, Real k_perp_sq) const override;
	Real xh(unsigned fl, Real x, Real Q_sq, Real k_perp_sq) const override;
	Real xhL(unsigned fl, Real x, Real Q_sq, Real k_perp_sq) const override;
	Real xhT(unsigned fl, Real x, Real Q_sq, Real k_perp_sq) const override;
	Real xhTperp(unsigned fl, Real x, Real Q_sq, Real k_perp_sq) const override;
	Real xe(unsigned fl, Real x, Real Q_sq, Real k_perp_sq) const override;
	Real xeL(unsigned fl, Real x, Real Q_sq, Real k_perp_sq) const override;
	Real xeT(unsigned fl, Real x, Real Q_sq, Real k_perp_sq) const override;
	Real xeTperp(unsigned fl, Real x, Real Q_sq, Real k_perp_sq) const override;

	Real D1(constant::Hadron h, unsigned fl, Real z, Real Q_sq, Real p_perp_sq) const override;
	Real H1perp(constant::Hadron h, unsigned fl, Real z, Real Q_sq, Real p_perp_sq) const override;
	Real Dperp_tilde(constant::Hadron h, unsigned fl, Real z, Real Q_sq, Real p_perp_sq) const override;
	Real H_tilde(constant::Hadron h, unsigned fl, Real z, Real Q_sq, Real p_perp_sq) const override;
	Real Gperp_tilde(constant::Hadron h, unsigned fl, Real z, Real Q_sq, Real p_perp_sq) const override;
	Real E_tilde(constant::Hadron h, unsigned fl, Real z, Real Q_sq, Real p_perp_sq) const override;
};

class TmdWwSet : public TmdSet {
public:
	TmdWwSet(
		unsigned flavor_count,
		constant::Nucleus target) :
		TmdSet(
			flavor_count,
			target) { }
	virtual ~TmdWwSet() = default;

	virtual Real xf1(unsigned fl, Real x, Real Q_sq, Real k_perp_sq) const;
	virtual Real xf1Tperp(unsigned fl, Real x, Real Q_sq, Real k_perp_sq) const;
	virtual Real xg1(unsigned fl, Real x, Real Q_sq, Real k_perp_sq) const;
	virtual Real xg1Tperp(unsigned fl, Real x, Real Q_sq, Real k_perp_sq) const;
	virtual Real xh1(unsigned fl, Real x, Real Q_sq, Real k_perp_sq) const;
	virtual Real xh1perp(unsigned fl, Real x, Real Q_sq, Real k_perp_sq) const;
	virtual Real xh1Lperp(unsigned fl, Real x, Real Q_sq, Real k_perp_sq) const;
	virtual Real xh1Tperp(unsigned fl, Real x, Real Q_sq, Real k_perp_sq) const;

	virtual Real D1(constant::Hadron h, unsigned fl, Real z, Real Q_sq, Real p_perp_sq) const;
	virtual Real H1perp(constant::Hadron h, unsigned fl, Real z, Real Q_sq, Real p_perp_sq) const;

	Real xf1TperpM1(unsigned fl, Real x, Real Q_sq, Real k_perp_sq) const;
	Real xg1TperpM1(unsigned fl, Real x, Real Q_sq, Real k_perp_sq) const;
	Real xh1perpM1(unsigned fl, Real x, Real Q_sq, Real k_perp_sq) const;
	Real xh1LperpM1(unsigned fl, Real x, Real Q_sq, Real k_perp_sq) const;
	Real xh1TperpM1(unsigned fl, Real x, Real Q_sq, Real k_perp_sq) const;

	Real xfT(unsigned fl, Real x, Real Q_sq, Real k_perp_sq) const override;
	Real xfperp(unsigned fl, Real x, Real Q_sq, Real k_perp_sq) const override;
	Real xfLperp(unsigned fl, Real x, Real Q_sq, Real k_perp_sq) const override;
	Real xfTperp(unsigned fl, Real x, Real Q_sq, Real k_perp_sq) const override;
	Real xgT(unsigned fl, Real x, Real Q_sq, Real k_perp_sq) const override;
	Real xgperp(unsigned fl, Real x, Real Q_sq, Real k_perp_sq) const override;
	Real xgLperp(unsigned fl, Real x, Real Q_sq, Real k_perp_sq) const override;
	Real xgTperp(unsigned fl, Real x, Real Q_sq, Real k_perp_sq) const override;
	Real xh(unsigned fl, Real x, Real Q_sq, Real k_perp_sq) const override;
	Real xhL(unsigned fl, Real x, Real Q_sq, Real k_perp_sq) const override;
	Real xhT(unsigned fl, Real x, Real Q_sq, Real k_perp_sq) const override;
	Real xhTperp(unsigned fl, Real x, Real Q_sq, Real k_perp_sq) const override;
	Real xe(unsigned fl, Real x, Real Q_sq, Real k_perp_sq) const override;
	Real xeL(unsigned fl, Real x, Real Q_sq, Real k_perp_sq) const override;
	Real xeT(unsigned fl, Real x, Real Q_sq, Real k_perp_sq) const override;
	Real xeTperp(unsigned fl, Real x, Real Q_sq, Real k_perp_sq) const override;

	Real Dperp_tilde(constant::Hadron h, unsigned fl, Real z, Real Q_sq, Real p_perp_sq) const override;
	Real H_tilde(constant::Hadron h, unsigned fl, Real z, Real Q_sq, Real p_perp_sq) const override;
	Real Gperp_tilde(constant::Hadron h, unsigned fl, Real z, Real Q_sq, Real p_perp_sq) const override;
	Real E_tilde(constant::Hadron h, unsigned fl, Real z, Real Q_sq, Real p_perp_sq) const override;
};

class TmdGaussianWwSet : public TmdGaussianSet {
public:
	TmdGaussianWwSet(
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
		Real mean_H1perp);
	virtual ~TmdGaussianWwSet() = default;

	virtual Real xf1(unsigned fl, Real x, Real Q_sq) const;
	virtual Real xf1Tperp(unsigned fl, Real x, Real Q_sq) const;
	virtual Real xg1(unsigned fl, Real x, Real Q_sq) const;
	virtual Real xg1Tperp(unsigned fl, Real x, Real Q_sq) const;
	virtual Real xh1(unsigned fl, Real x, Real Q_sq) const;
	virtual Real xh1perp(unsigned fl, Real x, Real Q_sq) const;
	virtual Real xh1Lperp(unsigned fl, Real x, Real Q_sq) const;
	virtual Real xh1Tperp(unsigned fl, Real x, Real Q_sq) const;

	virtual Real D1(constant::Hadron h, unsigned fl, Real z, Real Q_sq) const;
	virtual Real H1perp(constant::Hadron h, unsigned fl, Real z, Real Q_sq) const;

	Real xf1TperpM1(unsigned fl, Real x, Real Q_sq) const;
	Real xg1TperpM1(unsigned fl, Real x, Real Q_sq) const;
	Real xh1perpM1(unsigned fl, Real x, Real Q_sq) const;
	Real xh1LperpM1(unsigned fl, Real x, Real Q_sq) const;
	Real xh1TperpM1(unsigned fl, Real x, Real Q_sq) const;

	Real xfT(unsigned fl, Real x, Real Q_sq) const override;
	Real xfperp(unsigned fl, Real x, Real Q_sq) const override;
	Real xfLperp(unsigned fl, Real x, Real Q_sq) const override;
	Real xfTperp(unsigned fl, Real x, Real Q_sq) const override;
	Real xgT(unsigned fl, Real x, Real Q_sq) const override;
	Real xgperp(unsigned fl, Real x, Real Q_sq) const override;
	Real xgLperp(unsigned fl, Real x, Real Q_sq) const override;
	Real xgTperp(unsigned fl, Real x, Real Q_sq) const override;
	Real xh(unsigned fl, Real x, Real Q_sq) const override;
	Real xhL(unsigned fl, Real x, Real Q_sq) const override;
	Real xhT(unsigned fl, Real x, Real Q_sq) const override;
	Real xhTperp(unsigned fl, Real x, Real Q_sq) const override;
	Real xe(unsigned fl, Real x, Real Q_sq) const override;
	Real xeL(unsigned fl, Real x, Real Q_sq) const override;
	Real xeT(unsigned fl, Real x, Real Q_sq) const override;
	Real xeTperp(unsigned fl, Real x, Real Q_sq) const override;

	Real Dperp_tilde(constant::Hadron h, unsigned fl, Real z, Real Q_sq) const override;
	Real H_tilde(constant::Hadron h, unsigned fl, Real z, Real Q_sq) const override;
	Real Gperp_tilde(constant::Hadron h, unsigned fl, Real z, Real Q_sq) const override;
	Real E_tilde(constant::Hadron h, unsigned fl, Real z, Real Q_sq) const override;
};

}
}

#endif

