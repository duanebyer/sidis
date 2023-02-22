#ifndef SIDIS_TMD_HPP
#define SIDIS_TMD_HPP

#include <vector>

#include "sidis/particle.hpp"

namespace sidis {
namespace sf {

/**
 * \defgroup TmdGroup Transverse momentum distributions (TMD)
 * Interfaces for user-defined transverse momentum distributions (TMDs).
 */
/// \{

/**
 * Complete set of transverse momentum distributions (TMDs) and also
 * fragmentation functions (FFs) bundled together. This abstract class is to be
 * derived for user-provided TMDs and FFs. The TMDs and FFs are indexed by a
 * flavor index, where the number of supported flavors is also user-provided.
 *
 * When possible, it is recommended to derive GaussianTmdSet, WwTmdSet, or
 * GaussianWwTmdSet instead.
 */
class TmdSet {
public:
	/// The number of flavors supported by the TmdSet.
	unsigned const flavor_count;
	/// What type of part::Nucleus the TMDs are valid for.
	part::Nucleus const target;

	/// Initialize a TMDSet with \p flavor_count number of flavors and for the
	/// specified target.
	TmdSet(unsigned flavor_count, part::Nucleus target) :
		flavor_count(flavor_count),
		target(target) { }
	TmdSet(TmdSet const&) = delete;
	TmdSet(TmdSet&&) = delete;
	TmdSet& operator=(TmdSet const&) = delete;
	TmdSet& operator=(TmdSet&&) = delete;
	virtual ~TmdSet() = default;

	/// Charge of a given flavor (for use as the weighting when computing
	/// structure functions).
	virtual Real charge(unsigned fl) const = 0;

	/// \name Transverse momentum distributions
	/// By default, these return zero.
	/// \{
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
	/// \}

	/// \name Fragmentation functions
	/// By default, these return zero.
	/// \{
	virtual Real D1(part::Hadron h, unsigned fl, Real z, Real Q_sq, Real p_perp_sq) const;
	virtual Real H1perp(part::Hadron h, unsigned fl, Real z, Real Q_sq, Real p_perp_sq) const;
	virtual Real Dperp_tilde(part::Hadron h, unsigned fl, Real z, Real Q_sq, Real p_perp_sq) const;
	virtual Real H_tilde(part::Hadron h, unsigned fl, Real z, Real Q_sq, Real p_perp_sq) const;
	virtual Real Gperp_tilde(part::Hadron h, unsigned fl, Real z, Real Q_sq, Real p_perp_sq) const;
	virtual Real E_tilde(part::Hadron h, unsigned fl, Real z, Real Q_sq, Real p_perp_sq) const;
	/// \}
};

/**
 * Gaussian TMDs and FFs for more efficient structure function computation.
 *
 * \sa TmdSet
 */
class GaussianTmdSet : public TmdSet {
public:
	/// \name Gaussian widths of TMDs
	/// \{
	std::vector<Real> const mean_f1;
	std::vector<Real> const mean_f1Tperp;
	std::vector<Real> const mean_fT;
	std::vector<Real> const mean_fperp;
	std::vector<Real> const mean_fLperp;
	std::vector<Real> const mean_fTperp;
	std::vector<Real> const mean_g1;
	std::vector<Real> const mean_g1Tperp;
	std::vector<Real> const mean_gT;
	std::vector<Real> const mean_gperp;
	std::vector<Real> const mean_gLperp;
	std::vector<Real> const mean_gTperp;
	std::vector<Real> const mean_h1;
	std::vector<Real> const mean_h1perp;
	std::vector<Real> const mean_h1Lperp;
	std::vector<Real> const mean_h1Tperp;
	std::vector<Real> const mean_h;
	std::vector<Real> const mean_hL;
	std::vector<Real> const mean_hT;
	std::vector<Real> const mean_hTperp;
	std::vector<Real> const mean_e;
	std::vector<Real> const mean_eL;
	std::vector<Real> const mean_eT;
	std::vector<Real> const mean_eTperp;
	/// \}

	/// \name Gaussian widths of FFs
	/// \{
	std::vector<Real> const mean_D1;
	std::vector<Real> const mean_H1perp;
	std::vector<Real> const mean_Dperp_tilde;
	std::vector<Real> const mean_H_tilde;
	std::vector<Real> const mean_Gperp_tilde;
	std::vector<Real> const mean_E_tilde;
	/// \}

	/// Initialize a GaussianTmdSet with the provided Gaussian widths.
	GaussianTmdSet(
		unsigned flavor_count,
		part::Nucleus target,
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
	/// Initialize a GaussianTmdSet with the provided flavor-dependent widths.
	GaussianTmdSet(
		unsigned flavor_count,
		part::Nucleus target,
		Real const* mean_f1,
		Real const* mean_f1Tperp,
		Real const* mean_fT,
		Real const* mean_fperp,
		Real const* mean_fLperp,
		Real const* mean_fTperp,
		Real const* mean_g1,
		Real const* mean_g1Tperp,
		Real const* mean_gT,
		Real const* mean_gperp,
		Real const* mean_gLperp,
		Real const* mean_gTperp,
		Real const* mean_h1,
		Real const* mean_h1perp,
		Real const* mean_h1Lperp,
		Real const* mean_h1Tperp,
		Real const* mean_h,
		Real const* mean_hL,
		Real const* mean_hT,
		Real const* mean_hTperp,
		Real const* mean_e,
		Real const* mean_eL,
		Real const* mean_eT,
		Real const* mean_eTperp,
		Real const* mean_D1,
		Real const* mean_H1perp,
		Real const* mean_Dperp_tilde,
		Real const* mean_H_tilde,
		Real const* mean_Gperp_tilde,
		Real const* mean_E_tilde);
	virtual ~GaussianTmdSet() = default;

	/// \name Reduced transverse momentum distributions
	/// These TMDs do not have explicit \f$\pmb{k}_{\perp}\f$ dependence, since
	/// the \f$\pmb{k}_{\perp}\f$ dependence is given by the Gaussian widths.
	/// \{
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
	/// \}

	/// \name Reduced fragmentation functions
	/// These FFs do not have explicit \f$\pmb{P}_{\perp}\f$ dependence, since
	/// the \f$\pmb{P}_{\perp}\f$ dependence is given by the Gaussian widths.
	/// \{
	virtual Real D1(part::Hadron h, unsigned fl, Real z, Real Q_sq) const;
	virtual Real H1perp(part::Hadron h, unsigned fl, Real z, Real Q_sq) const;
	virtual Real Dperp_tilde(part::Hadron h, unsigned fl, Real z, Real Q_sq) const;
	virtual Real H_tilde(part::Hadron h, unsigned fl, Real z, Real Q_sq) const;
	virtual Real Gperp_tilde(part::Hadron h, unsigned fl, Real z, Real Q_sq) const;
	virtual Real E_tilde(part::Hadron h, unsigned fl, Real z, Real Q_sq) const;
	/// \}

	/// \name Transverse momentum distributions
	/// \{
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
	/// \}

	/// \name Fragmentation functions
	/// \{
	Real D1(part::Hadron h, unsigned fl, Real z, Real Q_sq, Real p_perp_sq) const override;
	Real H1perp(part::Hadron h, unsigned fl, Real z, Real Q_sq, Real p_perp_sq) const override;
	Real Dperp_tilde(part::Hadron h, unsigned fl, Real z, Real Q_sq, Real p_perp_sq) const override;
	Real H_tilde(part::Hadron h, unsigned fl, Real z, Real Q_sq, Real p_perp_sq) const override;
	Real Gperp_tilde(part::Hadron h, unsigned fl, Real z, Real Q_sq, Real p_perp_sq) const override;
	Real E_tilde(part::Hadron h, unsigned fl, Real z, Real Q_sq, Real p_perp_sq) const override;
	/// \}
};

/**
 * Use Wandzura-Wilczek-type (WW-type) approximations compute all TMDs and FFs
 * in terms of a small basis of TMDs and FFs. See \cite bastami2019ww.
 *
 * \sa TmdSet
 */
class WwTmdSet : public TmdSet {
public:
	/// Initialize a WwTMDSet with \p flavor_count number of flavors and for the
	/// specified target.
	WwTmdSet(
		unsigned flavor_count,
		part::Nucleus target) :
		TmdSet(
			flavor_count,
			target) { }
	virtual ~WwTmdSet() = default;

	/// \name Base set of transverse momentum distributions
	/// By default, these return zero.
	/// \{
	virtual Real xf1(unsigned fl, Real x, Real Q_sq, Real k_perp_sq) const override;
	virtual Real xf1Tperp(unsigned fl, Real x, Real Q_sq, Real k_perp_sq) const override;
	virtual Real xg1(unsigned fl, Real x, Real Q_sq, Real k_perp_sq) const override;
	virtual Real xg1Tperp(unsigned fl, Real x, Real Q_sq, Real k_perp_sq) const override;
	virtual Real xh1(unsigned fl, Real x, Real Q_sq, Real k_perp_sq) const override;
	virtual Real xh1perp(unsigned fl, Real x, Real Q_sq, Real k_perp_sq) const override;
	virtual Real xh1Lperp(unsigned fl, Real x, Real Q_sq, Real k_perp_sq) const override;
	virtual Real xh1Tperp(unsigned fl, Real x, Real Q_sq, Real k_perp_sq) const override;
	/// \}

	/// \name Base set of fragmentation functions
	/// By default, these return zero.
	/// \{
	virtual Real D1(part::Hadron h, unsigned fl, Real z, Real Q_sq, Real p_perp_sq) const override;
	virtual Real H1perp(part::Hadron h, unsigned fl, Real z, Real Q_sq, Real p_perp_sq) const override;
	/// \}

	/// \name Transverse momentum distribution moments
	/// Important first moments of TMDs.
	/// \{
	Real xf1TperpM1(unsigned fl, Real x, Real Q_sq, Real k_perp_sq) const;
	Real xg1TperpM1(unsigned fl, Real x, Real Q_sq, Real k_perp_sq) const;
	Real xh1perpM1(unsigned fl, Real x, Real Q_sq, Real k_perp_sq) const;
	Real xh1LperpM1(unsigned fl, Real x, Real Q_sq, Real k_perp_sq) const;
	Real xh1TperpM1(unsigned fl, Real x, Real Q_sq, Real k_perp_sq) const;
	/// \}

	/// \name Transverse momentum distributions
	/// \{
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
	/// \}

	/// \name Fragmentation functions
	/// \{
	Real Dperp_tilde(part::Hadron h, unsigned fl, Real z, Real Q_sq, Real p_perp_sq) const override;
	Real H_tilde(part::Hadron h, unsigned fl, Real z, Real Q_sq, Real p_perp_sq) const override;
	Real Gperp_tilde(part::Hadron h, unsigned fl, Real z, Real Q_sq, Real p_perp_sq) const override;
	Real E_tilde(part::Hadron h, unsigned fl, Real z, Real Q_sq, Real p_perp_sq) const override;
	/// \}
};

/**
 * Combine Gaussian and WW-type approximations together. Because of some
 * inconsistencies between Gaussian and WW-type approximations, the method of
 * combining these two (take from \cite bastami2019ww) is to first apply the
 * Gaussian approximation when performing convolutions, and then to apply the
 * WW-type approximation to the reduced TMDs and FFs without
 * \f$\pmb{k}_{\perp}\f$ or \f$\pmb{P}_{\perp}\f$ dependence.
 *
 * \sa TmdSet
 */
class GaussianWwTmdSet : public GaussianTmdSet {
public:
	/// Initialize a GaussianWwTmdSet with the provided Gaussian widths.
	GaussianWwTmdSet(
		unsigned flavor_count,
		part::Nucleus target,
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
	/// Initialize a GaussianWwTmdSet with the provided flavor-dependent widths.
	GaussianWwTmdSet(
		unsigned flavor_count,
		part::Nucleus target,
		Real const* mean_f1,
		Real const* mean_f1Tperp,
		Real const* mean_fT,
		Real const* mean_g1,
		Real const* mean_g1Tperp,
		Real const* mean_gT,
		Real const* mean_h1,
		Real const* mean_h1perp,
		Real const* mean_h1Lperp,
		Real const* mean_h1Tperp,
		Real const* mean_h,
		Real const* mean_hL,
		Real const* mean_hT,
		Real const* mean_hTperp,
		Real const* mean_D1,
		Real const* mean_H1perp);
	virtual ~GaussianWwTmdSet() = default;

	/// \name Base set of reduced transverse momentum distributions
	/// These TMDs do not have explicit \f$\pmb{k}_{\perp}\f$ dependence, since
	/// the \f$\pmb{k}_{\perp}\f$ dependence is given by the Gaussian widths.
	/// \{
	virtual Real xf1(unsigned fl, Real x, Real Q_sq) const override;
	virtual Real xf1Tperp(unsigned fl, Real x, Real Q_sq) const override;
	virtual Real xg1(unsigned fl, Real x, Real Q_sq) const override;
	virtual Real xg1Tperp(unsigned fl, Real x, Real Q_sq) const override;
	virtual Real xh1(unsigned fl, Real x, Real Q_sq) const override;
	virtual Real xh1perp(unsigned fl, Real x, Real Q_sq) const override;
	virtual Real xh1Lperp(unsigned fl, Real x, Real Q_sq) const override;
	virtual Real xh1Tperp(unsigned fl, Real x, Real Q_sq) const override;
	/// \}

	/// \name Base set of reduced fragmentation functions
	/// These FFs do not have explicit \f$\pmb{P}_{\perp}\f$ dependence, since
	/// the \f$\pmb{P}_{\perp}\f$ dependence is given by the Gaussian widths.
	/// \{
	virtual Real D1(part::Hadron h, unsigned fl, Real z, Real Q_sq) const override;
	virtual Real H1perp(part::Hadron h, unsigned fl, Real z, Real Q_sq) const override;
	/// \}

	/// \name Transverse momentum distribution moments
	/// Important reduced first moments of TMDs.
	/// \{
	Real xf1TperpM1(unsigned fl, Real x, Real Q_sq) const;
	Real xg1TperpM1(unsigned fl, Real x, Real Q_sq) const;
	Real xh1perpM1(unsigned fl, Real x, Real Q_sq) const;
	Real xh1LperpM1(unsigned fl, Real x, Real Q_sq) const;
	Real xh1TperpM1(unsigned fl, Real x, Real Q_sq) const;

	/// Transverse momentum distributions
	/// \{
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
	/// \}

	/// Fragmentation functions
	/// \{
	Real Dperp_tilde(part::Hadron h, unsigned fl, Real z, Real Q_sq) const override;
	Real H_tilde(part::Hadron h, unsigned fl, Real z, Real Q_sq) const override;
	Real Gperp_tilde(part::Hadron h, unsigned fl, Real z, Real Q_sq) const override;
	Real E_tilde(part::Hadron h, unsigned fl, Real z, Real Q_sq) const override;
	/// \}
};
/// \}

}
}

#endif

