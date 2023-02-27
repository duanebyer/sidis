#ifndef SIDIS_TMD_HPP
#define SIDIS_TMD_HPP

#include "sidis/flavor_vec.hpp"
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
	/// The charges of each flavor.
	FlavorVec const charges;
	/// What type of part::Nucleus the TMDs are valid for.
	part::Nucleus const target;

	/// Initialize a TMDSet with \p flavor_count number of flavors and for the
	/// specified target.
	TmdSet(unsigned flavor_count, FlavorVec const& charges, part::Nucleus target);
	TmdSet(TmdSet const&) = delete;
	TmdSet(TmdSet&&) = delete;
	TmdSet& operator=(TmdSet const&) = delete;
	TmdSet& operator=(TmdSet&&) = delete;
	virtual ~TmdSet() = default;

	/// \name Transverse momentum distributions
	/// By default, these return zero.
	/// \{
	virtual FlavorVec xf1(Real x, Real Q_sq, Real k_perp_sq) const;
	virtual FlavorVec xf1Tperp(Real x, Real Q_sq, Real k_perp_sq) const;
	virtual FlavorVec xfT(Real x, Real Q_sq, Real k_perp_sq) const;
	virtual FlavorVec xfperp(Real x, Real Q_sq, Real k_perp_sq) const;
	virtual FlavorVec xfLperp(Real x, Real Q_sq, Real k_perp_sq) const;
	virtual FlavorVec xfTperp(Real x, Real Q_sq, Real k_perp_sq) const;
	virtual FlavorVec xg1(Real x, Real Q_sq, Real k_perp_sq) const;
	virtual FlavorVec xg1Tperp(Real x, Real Q_sq, Real k_perp_sq) const;
	virtual FlavorVec xgT(Real x, Real Q_sq, Real k_perp_sq) const;
	virtual FlavorVec xgperp(Real x, Real Q_sq, Real k_perp_sq) const;
	virtual FlavorVec xgLperp(Real x, Real Q_sq, Real k_perp_sq) const;
	virtual FlavorVec xgTperp(Real x, Real Q_sq, Real k_perp_sq) const;
	virtual FlavorVec xh1(Real x, Real Q_sq, Real k_perp_sq) const;
	virtual FlavorVec xh1perp(Real x, Real Q_sq, Real k_perp_sq) const;
	virtual FlavorVec xh1Lperp(Real x, Real Q_sq, Real k_perp_sq) const;
	virtual FlavorVec xh1Tperp(Real x, Real Q_sq, Real k_perp_sq) const;
	virtual FlavorVec xh(Real x, Real Q_sq, Real k_perp_sq) const;
	virtual FlavorVec xhL(Real x, Real Q_sq, Real k_perp_sq) const;
	virtual FlavorVec xhT(Real x, Real Q_sq, Real k_perp_sq) const;
	virtual FlavorVec xhTperp(Real x, Real Q_sq, Real k_perp_sq) const;
	virtual FlavorVec xe(Real x, Real Q_sq, Real k_perp_sq) const;
	virtual FlavorVec xeL(Real x, Real Q_sq, Real k_perp_sq) const;
	virtual FlavorVec xeT(Real x, Real Q_sq, Real k_perp_sq) const;
	virtual FlavorVec xeTperp(Real x, Real Q_sq, Real k_perp_sq) const;
	/// \}

	/// \name Fragmentation functions
	/// By default, these return zero.
	/// \{
	virtual FlavorVec D1(part::Hadron h, Real z, Real Q_sq, Real p_perp_sq) const;          // XX
	virtual FlavorVec H1perp(part::Hadron h, Real z, Real Q_sq, Real p_perp_sq) const;      // XX
	virtual FlavorVec Dperp_tilde(part::Hadron h, Real z, Real Q_sq, Real p_perp_sq) const; // UU UT LL LT
	virtual FlavorVec H_tilde(part::Hadron h, Real z, Real Q_sq, Real p_perp_sq) const;     // UX
	virtual FlavorVec Gperp_tilde(part::Hadron h, Real z, Real Q_sq, Real p_perp_sq) const; // UL UT LU LT
	virtual FlavorVec E_tilde(part::Hadron h, Real z, Real Q_sq, Real p_perp_sq) const;     // LX
	/// \}
};

struct GaussianWwTmdVars;

/**
 * Stores the width of the Gaussian TMDs and FFs, for use in GaussianTmdSet.
 * Technically, the square of the width (equivalent to the variance) is stored.
 * The variances can be be specified by flavor.
 *
 * \sa GaussianTmdSet
 * \sa GaussianWwTmdVars
 */
struct GaussianTmdVars {
	/// \name Variances of TMDs
	/// \{
	FlavorVec f1;
	FlavorVec f1Tperp;
	FlavorVec fT;
	FlavorVec fperp;
	FlavorVec fLperp;
	FlavorVec fTperp;
	FlavorVec g1;
	FlavorVec g1Tperp;
	FlavorVec gT;
	FlavorVec gperp;
	FlavorVec gLperp;
	FlavorVec gTperp;
	FlavorVec h1;
	FlavorVec h1perp;
	FlavorVec h1Lperp;
	FlavorVec h1Tperp;
	FlavorVec h;
	FlavorVec hL;
	FlavorVec hT;
	FlavorVec hTperp;
	FlavorVec e;
	FlavorVec eL;
	FlavorVec eT;
	FlavorVec eTperp;
	/// \}

	/// \name Variances of FFs
	/// \{
	FlavorVec D1;
	FlavorVec H1perp;
	FlavorVec Dperp_tilde;
	FlavorVec H_tilde;
	FlavorVec Gperp_tilde;
	FlavorVec E_tilde;
	/// \}

	/// Initializes the variances. Each variance is set to be infinity by
	/// default, which means the corresponding TMD will be completely neglected.
	explicit GaussianTmdVars();
	/// Initializes the variances using the variances for the Wandzura-Wilczek
	/// base set of TMDs and FFs.
	explicit GaussianTmdVars(GaussianWwTmdVars const& ww_vars);
};

/**
 * Gaussian TMDs and FFs for more efficient structure function computation.
 *
 * \sa TmdSet
 */
class GaussianTmdSet : public TmdSet {
public:
	/// Gaussian variances for each TMD.
	GaussianTmdVars const vars;

	/// Initialize a GaussianTmdSet with the provided Gaussian variances
	/// \p vars.
	GaussianTmdSet(
		unsigned flavor_count,
		FlavorVec const& charges,
		part::Nucleus target,
		GaussianTmdVars const& vars);
	virtual ~GaussianTmdSet() = default;

	/// \name Reduced transverse momentum distributions
	/// These TMDs do not have explicit \f$\pmb{k}_{\perp}\f$ dependence, since
	/// the \f$\pmb{k}_{\perp}\f$ dependence is given by the Gaussian widths.
	/// \{
	virtual FlavorVec xf1(Real x, Real Q_sq) const;
	virtual FlavorVec xf1Tperp(Real x, Real Q_sq) const;
	virtual FlavorVec xfT(Real x, Real Q_sq) const;
	virtual FlavorVec xfperp(Real x, Real Q_sq) const;
	virtual FlavorVec xfLperp(Real x, Real Q_sq) const;
	virtual FlavorVec xfTperp(Real x, Real Q_sq) const;
	virtual FlavorVec xg1(Real x, Real Q_sq) const;
	virtual FlavorVec xg1Tperp(Real x, Real Q_sq) const;
	virtual FlavorVec xgT(Real x, Real Q_sq) const;
	virtual FlavorVec xgperp(Real x, Real Q_sq) const;
	virtual FlavorVec xgLperp(Real x, Real Q_sq) const;
	virtual FlavorVec xgTperp(Real x, Real Q_sq) const;
	virtual FlavorVec xh1(Real x, Real Q_sq) const;
	virtual FlavorVec xh1perp(Real x, Real Q_sq) const;
	virtual FlavorVec xh1Lperp(Real x, Real Q_sq) const;
	virtual FlavorVec xh1Tperp(Real x, Real Q_sq) const;
	virtual FlavorVec xh(Real x, Real Q_sq) const;
	virtual FlavorVec xhL(Real x, Real Q_sq) const;
	virtual FlavorVec xhT(Real x, Real Q_sq) const;
	virtual FlavorVec xhTperp(Real x, Real Q_sq) const;
	virtual FlavorVec xe(Real x, Real Q_sq) const;
	virtual FlavorVec xeL(Real x, Real Q_sq) const;
	virtual FlavorVec xeT(Real x, Real Q_sq) const;
	virtual FlavorVec xeTperp(Real x, Real Q_sq) const;
	/// \}

	/// \name Reduced fragmentation functions
	/// These FFs do not have explicit \f$\pmb{P}_{\perp}\f$ dependence, since
	/// the \f$\pmb{P}_{\perp}\f$ dependence is given by the Gaussian widths.
	/// \{
	virtual FlavorVec D1(part::Hadron h, Real z, Real Q_sq) const;
	virtual FlavorVec H1perp(part::Hadron h, Real z, Real Q_sq) const;
	virtual FlavorVec Dperp_tilde(part::Hadron h, Real z, Real Q_sq) const;
	virtual FlavorVec H_tilde(part::Hadron h, Real z, Real Q_sq) const;
	virtual FlavorVec Gperp_tilde(part::Hadron h, Real z, Real Q_sq) const;
	virtual FlavorVec E_tilde(part::Hadron h, Real z, Real Q_sq) const;
	/// \}

	/// \name Transverse momentum distributions
	/// \{
	FlavorVec xf1(Real x, Real Q_sq, Real k_perp_sq) const override;
	FlavorVec xf1Tperp(Real x, Real Q_sq, Real k_perp_sq) const override;
	FlavorVec xfT(Real x, Real Q_sq, Real k_perp_sq) const override;
	FlavorVec xfperp(Real x, Real Q_sq, Real k_perp_sq) const override;
	FlavorVec xfLperp(Real x, Real Q_sq, Real k_perp_sq) const override;
	FlavorVec xfTperp(Real x, Real Q_sq, Real k_perp_sq) const override;
	FlavorVec xg1(Real x, Real Q_sq, Real k_perp_sq) const override;
	FlavorVec xg1Tperp(Real x, Real Q_sq, Real k_perp_sq) const override;
	FlavorVec xgT(Real x, Real Q_sq, Real k_perp_sq) const override;
	FlavorVec xgperp(Real x, Real Q_sq, Real k_perp_sq) const override;
	FlavorVec xgLperp(Real x, Real Q_sq, Real k_perp_sq) const override;
	FlavorVec xgTperp(Real x, Real Q_sq, Real k_perp_sq) const override;
	FlavorVec xh1(Real x, Real Q_sq, Real k_perp_sq) const override;
	FlavorVec xh1perp(Real x, Real Q_sq, Real k_perp_sq) const override;
	FlavorVec xh1Lperp(Real x, Real Q_sq, Real k_perp_sq) const override;
	FlavorVec xh1Tperp(Real x, Real Q_sq, Real k_perp_sq) const override;
	FlavorVec xh(Real x, Real Q_sq, Real k_perp_sq) const override;
	FlavorVec xhL(Real x, Real Q_sq, Real k_perp_sq) const override;
	FlavorVec xhT(Real x, Real Q_sq, Real k_perp_sq) const override;
	FlavorVec xhTperp(Real x, Real Q_sq, Real k_perp_sq) const override;
	FlavorVec xe(Real x, Real Q_sq, Real k_perp_sq) const override;
	FlavorVec xeL(Real x, Real Q_sq, Real k_perp_sq) const override;
	FlavorVec xeT(Real x, Real Q_sq, Real k_perp_sq) const override;
	FlavorVec xeTperp(Real x, Real Q_sq, Real k_perp_sq) const override;
	/// \}

	/// \name Fragmentation functions
	/// \{
	FlavorVec D1(part::Hadron h, Real z, Real Q_sq, Real p_perp_sq) const override;
	FlavorVec H1perp(part::Hadron h, Real z, Real Q_sq, Real p_perp_sq) const override;
	FlavorVec Dperp_tilde(part::Hadron h, Real z, Real Q_sq, Real p_perp_sq) const override;
	FlavorVec H_tilde(part::Hadron h, Real z, Real Q_sq, Real p_perp_sq) const override;
	FlavorVec Gperp_tilde(part::Hadron h, Real z, Real Q_sq, Real p_perp_sq) const override;
	FlavorVec E_tilde(part::Hadron h, Real z, Real Q_sq, Real p_perp_sq) const override;
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
		FlavorVec const& charges,
		part::Nucleus target) :
		TmdSet(flavor_count, charges, target) { }
	virtual ~WwTmdSet() = default;

	/// \name Base set of transverse momentum distributions
	/// By default, these return zero.
	/// \{
	virtual FlavorVec xf1(Real x, Real Q_sq, Real k_perp_sq) const override;
	virtual FlavorVec xf1Tperp(Real x, Real Q_sq, Real k_perp_sq) const override;
	virtual FlavorVec xg1(Real x, Real Q_sq, Real k_perp_sq) const override;
	virtual FlavorVec xg1Tperp(Real x, Real Q_sq, Real k_perp_sq) const override;
	virtual FlavorVec xh1(Real x, Real Q_sq, Real k_perp_sq) const override;
	virtual FlavorVec xh1perp(Real x, Real Q_sq, Real k_perp_sq) const override;
	virtual FlavorVec xh1Lperp(Real x, Real Q_sq, Real k_perp_sq) const override;
	virtual FlavorVec xh1Tperp(Real x, Real Q_sq, Real k_perp_sq) const override;
	/// \}

	/// \name Base set of fragmentation functions
	/// By default, these return zero.
	/// \{
	virtual FlavorVec D1(part::Hadron h, Real z, Real Q_sq, Real p_perp_sq) const override;
	virtual FlavorVec H1perp(part::Hadron h, Real z, Real Q_sq, Real p_perp_sq) const override;
	/// \}

	/// \name Transverse momentum distribution moments
	/// Important first moments of TMDs.
	/// \{
	FlavorVec xf1TperpM1(Real x, Real Q_sq, Real k_perp_sq) const;
	FlavorVec xg1TperpM1(Real x, Real Q_sq, Real k_perp_sq) const;
	FlavorVec xh1perpM1(Real x, Real Q_sq, Real k_perp_sq) const;
	FlavorVec xh1LperpM1(Real x, Real Q_sq, Real k_perp_sq) const;
	FlavorVec xh1TperpM1(Real x, Real Q_sq, Real k_perp_sq) const;
	/// \}

	/// \name Transverse momentum distributions
	/// \{
	FlavorVec xfT(Real x, Real Q_sq, Real k_perp_sq) const override;
	FlavorVec xfperp(Real x, Real Q_sq, Real k_perp_sq) const override;
	FlavorVec xfLperp(Real x, Real Q_sq, Real k_perp_sq) const override;
	FlavorVec xfTperp(Real x, Real Q_sq, Real k_perp_sq) const override;
	FlavorVec xgT(Real x, Real Q_sq, Real k_perp_sq) const override;
	FlavorVec xgperp(Real x, Real Q_sq, Real k_perp_sq) const override;
	FlavorVec xgLperp(Real x, Real Q_sq, Real k_perp_sq) const override;
	FlavorVec xgTperp(Real x, Real Q_sq, Real k_perp_sq) const override;
	FlavorVec xh(Real x, Real Q_sq, Real k_perp_sq) const override;
	FlavorVec xhL(Real x, Real Q_sq, Real k_perp_sq) const override;
	FlavorVec xhT(Real x, Real Q_sq, Real k_perp_sq) const override;
	FlavorVec xhTperp(Real x, Real Q_sq, Real k_perp_sq) const override;
	FlavorVec xe(Real x, Real Q_sq, Real k_perp_sq) const override;
	FlavorVec xeL(Real x, Real Q_sq, Real k_perp_sq) const override;
	FlavorVec xeT(Real x, Real Q_sq, Real k_perp_sq) const override;
	FlavorVec xeTperp(Real x, Real Q_sq, Real k_perp_sq) const override;
	/// \}

	/// \name Fragmentation functions
	/// \{
	FlavorVec Dperp_tilde(part::Hadron h, Real z, Real Q_sq, Real p_perp_sq) const override;
	FlavorVec H_tilde(part::Hadron h, Real z, Real Q_sq, Real p_perp_sq) const override;
	FlavorVec Gperp_tilde(part::Hadron h, Real z, Real Q_sq, Real p_perp_sq) const override;
	FlavorVec E_tilde(part::Hadron h, Real z, Real Q_sq, Real p_perp_sq) const override;
	/// \}
};

/**
 * Similar to GaussianTmdVars, but only for the subset of TMDs and FFs needed in
 * the Wandzura-Wilczek-type approximation.
 *
 * \sa GaussianWwTmdSet
 * \sa GaussianTmdVars
 */
struct GaussianWwTmdVars {
	/// \name Variances of TMDs
	/// \{
	FlavorVec f1;
	FlavorVec f1Tperp;
	FlavorVec fT;
	FlavorVec g1;
	FlavorVec g1Tperp;
	FlavorVec gT;
	FlavorVec h1;
	FlavorVec h1perp;
	FlavorVec h1Lperp;
	FlavorVec h1Tperp;
	FlavorVec h;
	FlavorVec hL;
	FlavorVec hT;
	FlavorVec hTperp;
	/// \}

	/// \name Variances of FFs
	/// \{
	FlavorVec D1;
	FlavorVec H1perp;
	/// \}

	/// Initializes the variances. Each variance is set to be infinity by
	/// default, which means the corresponding TMD will be completely neglected.
	explicit GaussianWwTmdVars();
};

/**
 * Combine Gaussian and WW-type approximations together. Because of some
 * inconsistencies between Gaussian and WW-type approximations, the method of
 * combining these two (taken from \cite bastami2019ww) is to first apply the
 * Gaussian approximation when performing convolutions, and then to apply the
 * WW-type approximation to the reduced TMDs and FFs without
 * \f$\pmb{k}_{\perp}\f$ or \f$\pmb{P}_{\perp}\f$ dependence.
 *
 * \sa TmdSet
 */
class GaussianWwTmdSet : public GaussianTmdSet {
public:
	/// Initialize a GaussianWwTmdSet with the provided Gaussian variances.
	GaussianWwTmdSet(
		unsigned flavor_count,
		FlavorVec const& charges,
		part::Nucleus target,
		GaussianWwTmdVars const& vars);
	virtual ~GaussianWwTmdSet() = default;

	/// \name Base set of reduced transverse momentum distributions
	/// These TMDs do not have explicit \f$\pmb{k}_{\perp}\f$ dependence, since
	/// the \f$\pmb{k}_{\perp}\f$ dependence is given by the Gaussian widths.
	/// \{
	virtual FlavorVec xf1(Real x, Real Q_sq) const override;
	virtual FlavorVec xf1Tperp(Real x, Real Q_sq) const override;
	virtual FlavorVec xg1(Real x, Real Q_sq) const override;
	virtual FlavorVec xg1Tperp(Real x, Real Q_sq) const override;
	virtual FlavorVec xh1(Real x, Real Q_sq) const override;
	virtual FlavorVec xh1perp(Real x, Real Q_sq) const override;
	virtual FlavorVec xh1Lperp(Real x, Real Q_sq) const override;
	virtual FlavorVec xh1Tperp(Real x, Real Q_sq) const override;
	/// \}

	/// \name Base set of reduced fragmentation functions
	/// These FFs do not have explicit \f$\pmb{P}_{\perp}\f$ dependence, since
	/// the \f$\pmb{P}_{\perp}\f$ dependence is given by the Gaussian widths.
	/// \{
	virtual FlavorVec D1(part::Hadron h, Real z, Real Q_sq) const override;
	virtual FlavorVec H1perp(part::Hadron h, Real z, Real Q_sq) const override;
	/// \}

	/// \name Transverse momentum distribution moments
	/// Important reduced first moments of TMDs.
	/// \{
	FlavorVec xf1TperpM1(Real x, Real Q_sq) const;
	FlavorVec xg1TperpM1(Real x, Real Q_sq) const;
	FlavorVec xh1perpM1(Real x, Real Q_sq) const;
	FlavorVec xh1LperpM1(Real x, Real Q_sq) const;
	FlavorVec xh1TperpM1(Real x, Real Q_sq) const;

	/// Transverse momentum distributions
	/// \{
	FlavorVec xfT(Real x, Real Q_sq) const override;
	FlavorVec xfperp(Real x, Real Q_sq) const override;
	FlavorVec xfLperp(Real x, Real Q_sq) const override;
	FlavorVec xfTperp(Real x, Real Q_sq) const override;
	FlavorVec xgT(Real x, Real Q_sq) const override;
	FlavorVec xgperp(Real x, Real Q_sq) const override;
	FlavorVec xgLperp(Real x, Real Q_sq) const override;
	FlavorVec xgTperp(Real x, Real Q_sq) const override;
	FlavorVec xh(Real x, Real Q_sq) const override;
	FlavorVec xhL(Real x, Real Q_sq) const override;
	FlavorVec xhT(Real x, Real Q_sq) const override;
	FlavorVec xhTperp(Real x, Real Q_sq) const override;
	FlavorVec xe(Real x, Real Q_sq) const override;
	FlavorVec xeL(Real x, Real Q_sq) const override;
	FlavorVec xeT(Real x, Real Q_sq) const override;
	FlavorVec xeTperp(Real x, Real Q_sq) const override;
	/// \}

	/// Fragmentation functions
	/// \{
	FlavorVec Dperp_tilde(part::Hadron h, Real z, Real Q_sq) const override;
	FlavorVec H_tilde(part::Hadron h, Real z, Real Q_sq) const override;
	FlavorVec Gperp_tilde(part::Hadron h, Real z, Real Q_sq) const override;
	FlavorVec E_tilde(part::Hadron h, Real z, Real Q_sq) const override;
	/// \}
};
/// \}

}
}

#endif

