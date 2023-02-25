#ifndef SIDIS_STRUCTURE_FUNCTION_HPP
#define SIDIS_STRUCTURE_FUNCTION_HPP

#include "sidis/numeric.hpp"
#include "sidis/particle.hpp"
#include "sidis/tmd.hpp"

namespace sidis {
namespace sf {

struct SfBaseUU;
struct SfBaseUL;
struct SfBaseUT;
struct SfBaseUP;
struct SfBaseLU;
struct SfBaseLL;
struct SfBaseLT;
struct SfBaseLP;

struct SfUU;
struct SfUL;
struct SfUT;
struct SfUP;
struct SfLU;
struct SfLL;
struct SfLT;
struct SfLP;

/**
 * \defgroup SfGroup Structure functions
 * Classes related to computing the structure functions.
 */
/// \{

/**
 * Complete set of structure functions bundled together. This abstract class is
 * to be derived for user-provided structure functions. If the structure
 * functions can be factorized into transverse momentum distributions (TMDs),
 * then one of the derived classes TmdSfSet, GaussianTmdSfSet, WwTmdSfSet, or
 * GaussianWwTmdSfSet will be more suitable.
 *
 * This class contains all leading twist and sub-leading twist structure
 * functions.
 */
class SfSet {
public:
	/// What type of part::Nucleus the structure functions are valid for.
	part::Nucleus const target;

	/// Initialize a SfSet for the specified target.
	SfSet(part::Nucleus target) : target(target) { }
	SfSet(SfSet const&) = delete;
	SfSet(SfSet&&) = delete;
	SfSet& operator=(SfSet const&) = delete;
	SfSet& operator=(SfSet&&) = delete;
	virtual ~SfSet() = default;

	/// \name Structure functions
	/// All leading twist and sub-leading twist SIDIS structure functions. By
	/// default, each of these returns zero.
	/// \{
	virtual Real F_UUL(part::Hadron h, Real x, Real z, Real Q_sq, Real ph_t_sq) const;
	virtual Real F_UUT(part::Hadron h, Real x, Real z, Real Q_sq, Real ph_t_sq) const;
	virtual Real F_UU_cos_phih(part::Hadron h, Real x, Real z, Real Q_sq, Real ph_t_sq) const;
	virtual Real F_UU_cos_2phih(part::Hadron h, Real x, Real z, Real Q_sq, Real ph_t_sq) const;

	virtual Real F_UL_sin_phih(part::Hadron h, Real x, Real z, Real Q_sq, Real ph_t_sq) const;
	virtual Real F_UL_sin_2phih(part::Hadron h, Real x, Real z, Real Q_sq, Real ph_t_sq) const;

	virtual Real F_UTL_sin_phih_m_phis(part::Hadron h, Real x, Real z, Real Q_sq, Real ph_t_sq) const;
	virtual Real F_UTT_sin_phih_m_phis(part::Hadron h, Real x, Real z, Real Q_sq, Real ph_t_sq) const;
	virtual Real F_UT_sin_2phih_m_phis(part::Hadron h, Real x, Real z, Real Q_sq, Real ph_t_sq) const;
	virtual Real F_UT_sin_3phih_m_phis(part::Hadron h, Real x, Real z, Real Q_sq, Real ph_t_sq) const;
	virtual Real F_UT_sin_phis(part::Hadron h, Real x, Real z, Real Q_sq, Real ph_t_sq) const;
	virtual Real F_UT_sin_phih_p_phis(part::Hadron h, Real x, Real z, Real Q_sq, Real ph_t_sq) const;

	virtual Real F_LU_sin_phih(part::Hadron h, Real x, Real z, Real Q_sq, Real ph_t_sq) const;

	virtual Real F_LL(part::Hadron h, Real x, Real z, Real Q_sq, Real ph_t_sq) const;
	virtual Real F_LL_cos_phih(part::Hadron h, Real x, Real z, Real Q_sq, Real ph_t_sq) const;

	virtual Real F_LT_cos_phih_m_phis(part::Hadron h, Real x, Real z, Real Q_sq, Real ph_t_sq) const;
	virtual Real F_LT_cos_2phih_m_phis(part::Hadron h, Real x, Real z, Real Q_sq, Real ph_t_sq) const;
	virtual Real F_LT_cos_phis(part::Hadron h, Real x, Real z, Real Q_sq, Real ph_t_sq) const;
	/// \}

	/// \name Base structure function combinations
	/// Convenience methods for retrieving groups of structure functions based
	/// on beam and target polarizations. For example, SfSet::sf_base_ut() will
	/// retrieve all structure functions of the form \f$F_{UT}^{\cdot}\f$.
	///
	/// By default, this method packages together independent calls to the
	/// individual structure functions. The method is declared virtual in case a
	/// more efficient implementation is possible.
	/// \{
	virtual SfBaseUU sf_base_uu(part::Hadron h, Real x, Real z, Real Q_sq, Real ph_t_sq) const;
	virtual SfBaseUL sf_base_ul(part::Hadron h, Real x, Real z, Real Q_sq, Real ph_t_sq) const;
	virtual SfBaseUT sf_base_ut(part::Hadron h, Real x, Real z, Real Q_sq, Real ph_t_sq) const;
	virtual SfBaseUP sf_base_up(part::Hadron h, Real x, Real z, Real Q_sq, Real ph_t_sq) const;
	virtual SfBaseLU sf_base_lu(part::Hadron h, Real x, Real z, Real Q_sq, Real ph_t_sq) const;
	virtual SfBaseLL sf_base_ll(part::Hadron h, Real x, Real z, Real Q_sq, Real ph_t_sq) const;
	virtual SfBaseLT sf_base_lt(part::Hadron h, Real x, Real z, Real Q_sq, Real ph_t_sq) const;
	virtual SfBaseLP sf_base_lp(part::Hadron h, Real x, Real z, Real Q_sq, Real ph_t_sq) const;
	/// \}

	/// \name Structure function combinations
	/// Similar to the base structure function combinations, except that these
	/// methods will provide unpolarized variants as well. For example,
	/// SfSet::sf_ut() will retrieve structure functions of the form
	/// \f$F_{UT}^{\cdot}\f$ as well as \f$F_{UU}^{\cdot}\f$. The unpolarized
	/// variants are needed for the complete polarized cross-section.
	///
	/// By default, this method packages together independent calls to the
	/// individual structure functions. The method is declared virtual in case a
	/// more efficient implementation is possible.
	/// \{
	virtual SfUU sf_uu(part::Hadron h, Real x, Real z, Real Q_sq, Real ph_t_sq) const;
	virtual SfUL sf_ul(part::Hadron h, Real x, Real z, Real Q_sq, Real ph_t_sq) const;
	virtual SfUT sf_ut(part::Hadron h, Real x, Real z, Real Q_sq, Real ph_t_sq) const;
	virtual SfUP sf_up(part::Hadron h, Real x, Real z, Real Q_sq, Real ph_t_sq) const;
	virtual SfLU sf_lu(part::Hadron h, Real x, Real z, Real Q_sq, Real ph_t_sq) const;
	virtual SfLL sf_ll(part::Hadron h, Real x, Real z, Real Q_sq, Real ph_t_sq) const;
	virtual SfLT sf_lt(part::Hadron h, Real x, Real z, Real Q_sq, Real ph_t_sq) const;
	virtual SfLP sf_lp(part::Hadron h, Real x, Real z, Real Q_sq, Real ph_t_sq) const;
	/// \}
};

/**
 * Wrapper type for computing structure functions directly from the TMDs. The
 * structure functions are computed through convolution of the TMDs with the
 * fragmentation functions (FFs). If possible, prefer using a more specific
 * class such as GaussianTmdSfSet, WwTmdSfSet, or GaussianWwTmdSfSet to allow
 * the convolution to be computed more efficiently.
 *
 * The convolution used to compute the structure functions is:
 *
 * \f{equation}{
 *     F = x \sum_a e_a^2 \int d^2 \pmb{k}_{\perp} d^2 \pmb{P}_{\perp} \delta^{(2)}(z\pmb{k}_{\perp} + \pmb{P}_{\perp} - \pmb{p}_t)\omega f^a(x, |\pmb{k}_{\perp}|^2) D^a(x, |\pmb{D}_{\perp}|^2)
 * \f}
 *
 * Where \f$e_a\f$ is the charge of the parton flavor, \f$\omega\f$ is a
 * weighting factor, \f$f^a\f$ is a TMD, and \f$D^a\f$ is a FF.
 *
 * \sa TmdSet
 */
class TmdSfSet final : public SfSet {
public:
	/// The underlying TmdSet.
	TmdSet const& tmd_set;

	/// Initialize a TmdSfSet using the TMDs provided in \p tmd_set.
	TmdSfSet(TmdSet const& tmd_set) :
		SfSet(tmd_set.target),
		tmd_set(tmd_set) { };

	Real F_UUL(part::Hadron h, Real x, Real z, Real Q_sq, Real ph_t_sq) const override;
	Real F_UUT(part::Hadron h, Real x, Real z, Real Q_sq, Real ph_t_sq) const override;
	Real F_UU_cos_phih(part::Hadron h, Real x, Real z, Real Q_sq, Real ph_t_sq) const override;
	Real F_UU_cos_2phih(part::Hadron h, Real x, Real z, Real Q_sq, Real ph_t_sq) const override;

	Real F_UL_sin_phih(part::Hadron h, Real x, Real z, Real Q_sq, Real ph_t_sq) const override;
	Real F_UL_sin_2phih(part::Hadron h, Real x, Real z, Real Q_sq, Real ph_t_sq) const override;

	Real F_UTL_sin_phih_m_phis(part::Hadron h, Real x, Real z, Real Q_sq, Real ph_t_sq) const override;
	Real F_UTT_sin_phih_m_phis(part::Hadron h, Real x, Real z, Real Q_sq, Real ph_t_sq) const override;
	Real F_UT_sin_2phih_m_phis(part::Hadron h, Real x, Real z, Real Q_sq, Real ph_t_sq) const override;
	Real F_UT_sin_3phih_m_phis(part::Hadron h, Real x, Real z, Real Q_sq, Real ph_t_sq) const override;
	Real F_UT_sin_phis(part::Hadron h, Real x, Real z, Real Q_sq, Real ph_t_sq) const override;
	Real F_UT_sin_phih_p_phis(part::Hadron h, Real x, Real z, Real Q_sq, Real ph_t_sq) const override;

	Real F_LU_sin_phih(part::Hadron h, Real x, Real z, Real Q_sq, Real ph_t_sq) const override;

	Real F_LL(part::Hadron h, Real x, Real z, Real Q_sq, Real ph_t_sq) const override;
	Real F_LL_cos_phih(part::Hadron h, Real x, Real z, Real Q_sq, Real ph_t_sq) const override;

	Real F_LT_cos_phih_m_phis(part::Hadron h, Real x, Real z, Real Q_sq, Real ph_t_sq) const override;
	Real F_LT_cos_2phih_m_phis(part::Hadron h, Real x, Real z, Real Q_sq, Real ph_t_sq) const override;
	Real F_LT_cos_phis(part::Hadron h, Real x, Real z, Real Q_sq, Real ph_t_sq) const override;
};

/**
 * More efficient version of TmdSfSet specialized for Gaussian TMDs. The
 * convolution can be evaluated analytically.
 *
 * \sa GaussianTmdSet
 */
class GaussianTmdSfSet final : public SfSet {
public:
	/// The underlying GaussianTmdSet.
	GaussianTmdSet const& tmd_set;

	/// Initialize a GaussianTmdSfSet using the TMDs provided in \p tmd_set.
	GaussianTmdSfSet(GaussianTmdSet const& tmd_set) :
		SfSet(tmd_set.target),
		tmd_set(tmd_set) { };

	Real F_UUL(part::Hadron h, Real x, Real z, Real Q_sq, Real ph_t_sq) const override;
	Real F_UUT(part::Hadron h, Real x, Real z, Real Q_sq, Real ph_t_sq) const override;
	Real F_UU_cos_phih(part::Hadron h, Real x, Real z, Real Q_sq, Real ph_t_sq) const override;
	Real F_UU_cos_2phih(part::Hadron h, Real x, Real z, Real Q_sq, Real ph_t_sq) const override;

	Real F_UL_sin_phih(part::Hadron h, Real x, Real z, Real Q_sq, Real ph_t_sq) const override;
	Real F_UL_sin_2phih(part::Hadron h, Real x, Real z, Real Q_sq, Real ph_t_sq) const override;

	Real F_UTL_sin_phih_m_phis(part::Hadron h, Real x, Real z, Real Q_sq, Real ph_t_sq) const override;
	Real F_UTT_sin_phih_m_phis(part::Hadron h, Real x, Real z, Real Q_sq, Real ph_t_sq) const override;
	Real F_UT_sin_2phih_m_phis(part::Hadron h, Real x, Real z, Real Q_sq, Real ph_t_sq) const override;
	Real F_UT_sin_3phih_m_phis(part::Hadron h, Real x, Real z, Real Q_sq, Real ph_t_sq) const override;
	Real F_UT_sin_phis(part::Hadron h, Real x, Real z, Real Q_sq, Real ph_t_sq) const override;
	Real F_UT_sin_phih_p_phis(part::Hadron h, Real x, Real z, Real Q_sq, Real ph_t_sq) const override;

	Real F_LU_sin_phih(part::Hadron h, Real x, Real z, Real Q_sq, Real ph_t_sq) const override;

	Real F_LL(part::Hadron h, Real x, Real z, Real Q_sq, Real ph_t_sq) const override;
	Real F_LL_cos_phih(part::Hadron h, Real x, Real z, Real Q_sq, Real ph_t_sq) const override;

	Real F_LT_cos_phih_m_phis(part::Hadron h, Real x, Real z, Real Q_sq, Real ph_t_sq) const override;
	Real F_LT_cos_2phih_m_phis(part::Hadron h, Real x, Real z, Real Q_sq, Real ph_t_sq) const override;
	Real F_LT_cos_phis(part::Hadron h, Real x, Real z, Real Q_sq, Real ph_t_sq) const override;
};

/**
 * More efficient version of TmdSfSet specialized for TMDs with the
 * Wandzura-Wilczek-type (WW-type) approximation applied. The WW-type
 * approximation allows certain TMDs to be approximated in terms of others.
 *
 * \sa WwTmdSet
 */
class WwTmdSfSet final : public SfSet {
public:
	/// The underlying WwTmdSet.
	WwTmdSet const& tmd_set;

	/// Initialize a WwTmdSfSet using the TMDs provided in \p tmd_set.
	WwTmdSfSet(WwTmdSet const& tmd_set) :
		SfSet(tmd_set.target),
		tmd_set(tmd_set) { }

	Real F_UUL(part::Hadron h, Real x, Real z, Real Q_sq, Real ph_t_sq) const override;
	Real F_UUT(part::Hadron h, Real x, Real z, Real Q_sq, Real ph_t_sq) const override;
	Real F_UU_cos_phih(part::Hadron h, Real x, Real z, Real Q_sq, Real ph_t_sq) const override;
	Real F_UU_cos_2phih(part::Hadron h, Real x, Real z, Real Q_sq, Real ph_t_sq) const override;

	Real F_UL_sin_phih(part::Hadron h, Real x, Real z, Real Q_sq, Real ph_t_sq) const override;
	Real F_UL_sin_2phih(part::Hadron h, Real x, Real z, Real Q_sq, Real ph_t_sq) const override;

	Real F_UTL_sin_phih_m_phis(part::Hadron h, Real x, Real z, Real Q_sq, Real ph_t_sq) const override;
	Real F_UTT_sin_phih_m_phis(part::Hadron h, Real x, Real z, Real Q_sq, Real ph_t_sq) const override;
	Real F_UT_sin_2phih_m_phis(part::Hadron h, Real x, Real z, Real Q_sq, Real ph_t_sq) const override;
	Real F_UT_sin_3phih_m_phis(part::Hadron h, Real x, Real z, Real Q_sq, Real ph_t_sq) const override;
	Real F_UT_sin_phis(part::Hadron h, Real x, Real z, Real Q_sq, Real ph_t_sq) const override;
	Real F_UT_sin_phih_p_phis(part::Hadron h, Real x, Real z, Real Q_sq, Real ph_t_sq) const override;

	Real F_LU_sin_phih(part::Hadron h, Real x, Real z, Real Q_sq, Real ph_t_sq) const override;

	Real F_LL(part::Hadron h, Real x, Real z, Real Q_sq, Real ph_t_sq) const override;
	Real F_LL_cos_phih(part::Hadron h, Real x, Real z, Real Q_sq, Real ph_t_sq) const override;

	Real F_LT_cos_phih_m_phis(part::Hadron h, Real x, Real z, Real Q_sq, Real ph_t_sq) const override;
	Real F_LT_cos_2phih_m_phis(part::Hadron h, Real x, Real z, Real Q_sq, Real ph_t_sq) const override;
	Real F_LT_cos_phis(part::Hadron h, Real x, Real z, Real Q_sq, Real ph_t_sq) const override;
};

/**
 * More efficient version of TmdSfSet specialized for combining the Gaussian and
 * WW-type approximations. Since there are some issues when combining these two
 * approximations (see \cite bastami2019ww), this case must be treated
 * separately from GaussianTmdSfSet and WwTmdSfSet.
 *
 * \sa GaussianWwTmdSet
 */
class GaussianWwTmdSfSet final : public SfSet {
public:
	/// The underlying GaussianWwTmdSet.
	GaussianWwTmdSet const& tmd_set;

	/// Initialize a GaussianWwTmdSfSet using the TMDs provided in \p tmd_set.
	GaussianWwTmdSfSet(GaussianWwTmdSet const& tmd_set) :
		SfSet(tmd_set.target),
		tmd_set(tmd_set) { }

	Real F_UUL(part::Hadron h, Real x, Real z, Real Q_sq, Real ph_t_sq) const override;
	Real F_UUT(part::Hadron h, Real x, Real z, Real Q_sq, Real ph_t_sq) const override;
	Real F_UU_cos_phih(part::Hadron h, Real x, Real z, Real Q_sq, Real ph_t_sq) const override;
	Real F_UU_cos_2phih(part::Hadron h, Real x, Real z, Real Q_sq, Real ph_t_sq) const override;

	Real F_UL_sin_phih(part::Hadron h, Real x, Real z, Real Q_sq, Real ph_t_sq) const override;
	Real F_UL_sin_2phih(part::Hadron h, Real x, Real z, Real Q_sq, Real ph_t_sq) const override;

	Real F_UTL_sin_phih_m_phis(part::Hadron h, Real x, Real z, Real Q_sq, Real ph_t_sq) const override;
	Real F_UTT_sin_phih_m_phis(part::Hadron h, Real x, Real z, Real Q_sq, Real ph_t_sq) const override;
	Real F_UT_sin_2phih_m_phis(part::Hadron h, Real x, Real z, Real Q_sq, Real ph_t_sq) const override;
	Real F_UT_sin_3phih_m_phis(part::Hadron h, Real x, Real z, Real Q_sq, Real ph_t_sq) const override;
	Real F_UT_sin_phis(part::Hadron h, Real x, Real z, Real Q_sq, Real ph_t_sq) const override;
	Real F_UT_sin_phih_p_phis(part::Hadron h, Real x, Real z, Real Q_sq, Real ph_t_sq) const override;

	Real F_LU_sin_phih(part::Hadron h, Real x, Real z, Real Q_sq, Real ph_t_sq) const override;

	Real F_LL(part::Hadron h, Real x, Real z, Real Q_sq, Real ph_t_sq) const override;
	Real F_LL_cos_phih(part::Hadron h, Real x, Real z, Real Q_sq, Real ph_t_sq) const override;

	Real F_LT_cos_phih_m_phis(part::Hadron h, Real x, Real z, Real Q_sq, Real ph_t_sq) const override;
	Real F_LT_cos_2phih_m_phis(part::Hadron h, Real x, Real z, Real Q_sq, Real ph_t_sq) const override;
	Real F_LT_cos_phis(part::Hadron h, Real x, Real z, Real Q_sq, Real ph_t_sq) const override;
};

/**
 * \defgroup SfCombsBaseGroup Base structure function combinations
 * These types cache the results of structure function calculations from SfSet.
 * They are grouped by beam and target polarization. For instance, for
 * calculating the LT part of the cross-section through xs::born_base_lt(), you
 * would fill in the SfBaseLT type.
 *
 * To compute the full LT cross-section (including unpolarized contributions),
 * you would instead fill in the SfLT type.
 *
 * \sa SfSet
 */
/// \{
struct SfBaseUU {
	Real F_UUL;
	Real F_UUT;
	Real F_UU_cos_phih;
	Real F_UU_cos_2phih;
};
struct SfBaseUL {
	Real F_UL_sin_phih;
	Real F_UL_sin_2phih;
};
struct SfBaseUT {
	Real F_UTL_sin_phih_m_phis;
	Real F_UTT_sin_phih_m_phis;
	Real F_UT_sin_2phih_m_phis;
	Real F_UT_sin_3phih_m_phis;
	Real F_UT_sin_phis;
	Real F_UT_sin_phih_p_phis;
};
struct SfBaseUP {
	SfBaseUL ul;
	SfBaseUT ut;
};
struct SfBaseLU {
	Real F_LU_sin_phih;
};
struct SfBaseLL {
	Real F_LL;
	Real F_LL_cos_phih;
};
struct SfBaseLT {
	Real F_LT_cos_phih_m_phis;
	Real F_LT_cos_2phih_m_phis;
	Real F_LT_cos_phis;
};
struct SfBaseLP {
	SfBaseLL ll;
	SfBaseLT lt;
};
/// \}

/**
 * \defgroup SfCombsGroup Structure function combinations
 * These types combine the base structure function combinations together so that
 * a complete polarized cross-section can be calculated. For instance, if you
 * want to calculate an LT cross-section, you need to combine the UU, UT, LU,
 * and LT parts of the cross-section.
 *
 * \sa SfSet
 */
/// \{
struct SfUU {
	SfBaseUU uu;
};
struct SfUL {
	SfBaseUU uu;
	SfBaseUL ul;
};
struct SfUT {
	SfBaseUU uu;
	SfBaseUT ut;
};
struct SfUP {
	SfBaseUU uu;
	SfBaseUL ul;
	SfBaseUT ut;
};
struct SfLU {
	SfBaseUU uu;
	SfBaseLU lu;
};
struct SfLL {
	SfBaseUU uu;
	SfBaseUL ul;
	SfBaseLU lu;
	SfBaseLL ll;
};
struct SfLT {
	SfBaseUU uu;
	SfBaseUT ut;
	SfBaseLU lu;
	SfBaseLT lt;
};
struct SfLP {
	SfBaseUU uu;
	SfBaseUL ul;
	SfBaseUT ut;
	SfBaseLU lu;
	SfBaseLL ll;
	SfBaseLT lt;
};
/// \}
/// \}

}
}

#endif

