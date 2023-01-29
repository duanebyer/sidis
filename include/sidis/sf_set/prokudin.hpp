#ifndef SIDIS_SF_SET_PROKUDIN_HPP
#define SIDIS_SF_SET_PROKUDIN_HPP

#include "sidis/structure_function.hpp"
#include "sidis/tmd.hpp"

namespace sidis {

namespace constant {
	enum class Hadron;
}

namespace sf {
namespace set {

/**
 * Use data files from \cite bastami2019ww to calculate TMDs and FFs.
 */
class ProkudinTmdSet final : public GaussianWwTmdSet {
private:
	struct Impl;
	Impl* _impl;

public:
	ProkudinTmdSet();
	ProkudinTmdSet(ProkudinTmdSet const&) = delete;
	ProkudinTmdSet(ProkudinTmdSet&&) noexcept;
	ProkudinTmdSet& operator=(ProkudinTmdSet const&) = delete;
	ProkudinTmdSet& operator=(ProkudinTmdSet&& other) noexcept;
	virtual ~ProkudinTmdSet();

	Real charge(unsigned fl) const override;

	Real xf1(unsigned fl, Real x, Real Q_sq) const override;
	Real xf1Tperp(unsigned fl, Real x, Real Q_sq) const override;
	Real xg1(unsigned fl, Real x, Real Q_sq) const override;
	Real xg1Tperp(unsigned fl, Real x, Real Q_sq) const override;
	Real xh1(unsigned fl, Real x, Real Q_sq) const override;
	Real xh1perp(unsigned fl, Real x, Real Q_sq) const override;
	Real xh1Lperp(unsigned fl, Real x, Real Q_sq) const override;
	Real xh1Tperp(unsigned fl, Real x, Real Q_sq) const override;

	Real D1(part::Hadron h, unsigned fl, Real z, Real Q_sq) const override;
	Real H1perp(part::Hadron h, unsigned fl, Real z, Real Q_sq) const override;
};

/**
 * Use data files from \cite bastami2019ww to calculate structure functions.
 */
class ProkudinSfSet final : public SfSet {
	struct Impl;
	Impl* _impl;

	// Fragmentation functions.
	Real D1(part::Hadron h, unsigned fl, Real z, Real Q_sq) const;
	Real H1perpM1(part::Hadron h, unsigned fl, Real z, Real Q_sq) const;

	// Transverse momentum distributions.
	Real xf1(unsigned fl, Real x, Real Q_sq) const;
	Real xf1TperpM1(unsigned fl, Real x, Real Q_sq) const;
	Real xg1(unsigned fl, Real x, Real Q_sq) const;
	Real xgT(unsigned fl, Real x, Real Q_sq) const;
	Real xh1(unsigned fl, Real x, Real Q_sq) const;
	Real xh1M1(unsigned fl, Real x, Real Q_sq) const;
	Real xh1LperpM1(unsigned fl, Real x, Real Q_sq) const;
	Real xh1TperpM2(unsigned fl, Real x, Real Q_sq) const;
	Real xh1perpM1(unsigned fl, Real x, Real Q_sq) const;

public:
	ProkudinSfSet();
	ProkudinSfSet(ProkudinSfSet const& other) = delete;
	ProkudinSfSet(ProkudinSfSet&& other) noexcept;
	ProkudinSfSet& operator=(ProkudinSfSet const& other) = delete;
	ProkudinSfSet& operator=(ProkudinSfSet&& other) noexcept;
	virtual ~ProkudinSfSet();

	// Structure functions.
	Real F_UUT(part::Hadron h, Real x, Real z, Real Q_sq, Real ph_t_sq) const override;
	Real F_UU_cos_phih(part::Hadron h, Real x, Real z, Real Q_sq, Real ph_t_sq) const override;
	Real F_UU_cos_2phih(part::Hadron h, Real x, Real z, Real Q_sq, Real ph_t_sq) const override;

	Real F_UL_sin_phih(part::Hadron h, Real x, Real z, Real Q_sq, Real ph_t_sq) const override;
	Real F_UL_sin_2phih(part::Hadron h, Real x, Real z, Real Q_sq, Real ph_t_sq) const override;

	Real F_UTT_sin_phih_m_phis(part::Hadron h, Real x, Real z, Real Q_sq, Real ph_t_sq) const override;
	Real F_UT_sin_2phih_m_phis(part::Hadron h, Real x, Real z, Real Q_sq, Real ph_t_sq) const override;
	Real F_UT_sin_3phih_m_phis(part::Hadron h, Real x, Real z, Real Q_sq, Real ph_t_sq) const override;
	Real F_UT_sin_phis(part::Hadron h, Real x, Real z, Real Q_sq, Real ph_t_sq) const override;
	Real F_UT_sin_phih_p_phis(part::Hadron h, Real x, Real z, Real Q_sq, Real ph_t_sq) const override;

	Real F_LL(part::Hadron h, Real x, Real z, Real Q_sq, Real ph_t_sq) const override;
	Real F_LL_cos_phih(part::Hadron h, Real x, Real z, Real Q_sq, Real ph_t_sq) const override;

	Real F_LT_cos_phih_m_phis(part::Hadron h, Real x, Real z, Real Q_sq, Real ph_t_sq) const override;
	Real F_LT_cos_2phih_m_phis(part::Hadron h, Real x, Real z, Real Q_sq, Real ph_t_sq) const override;
	Real F_LT_cos_phis(part::Hadron h, Real x, Real z, Real Q_sq, Real ph_t_sq) const override;
};

}
}
}

#endif

