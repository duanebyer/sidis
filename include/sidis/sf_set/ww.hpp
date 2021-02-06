#ifndef SIDIS_SF_SET_WW_HPP
#define SIDIS_SF_SET_WW_HPP

#include <stdexcept>
#include <string>

#include "sidis/constant.hpp"
#include "sidis/numeric.hpp"
#include "sidis/structure_function.hpp"

namespace sidis {
namespace sf {
namespace set {

/**
 * Wandzura-Wilczek approximation to the structure functions [2].
 */
class WW final : public SfSet {
	struct Impl;
	Impl* _impl;

public:
	WW();
	~WW();

	WW(WW const& other) = delete;
	WW(WW&& other) noexcept;
	WW& operator=(WW const& other) = delete;
	WW& operator=(WW&& other) noexcept;

	// Structure functions.
	Real F_UUT(constant::Hadron h, Real x, Real z, Real Q_sq, Real ph_t_sq) const override;
	Real F_UU_cos_phih(constant::Hadron h, Real x, Real z, Real Q_sq, Real ph_t_sq) const override;
	Real F_UU_cos_2phih(constant::Hadron h, Real x, Real z, Real Q_sq, Real ph_t_sq) const override;

	Real F_UL_sin_phih(constant::Hadron h, Real x, Real z, Real Q_sq, Real ph_t_sq) const override;
	Real F_UL_sin_2phih(constant::Hadron h, Real x, Real z, Real Q_sq, Real ph_t_sq) const override;

	Real F_UTT_sin_phih_m_phis(constant::Hadron h, Real x, Real z, Real Q_sq, Real ph_t_sq) const override;
	Real F_UT_sin_2phih_m_phis(constant::Hadron h, Real x, Real z, Real Q_sq, Real ph_t_sq) const override;
	Real F_UT_sin_3phih_m_phis(constant::Hadron h, Real x, Real z, Real Q_sq, Real ph_t_sq) const override;
	Real F_UT_sin_phis(constant::Hadron h, Real x, Real z, Real Q_sq, Real ph_t_sq) const override;
	Real F_UT_sin_phih_p_phis(constant::Hadron h, Real x, Real z, Real Q_sq, Real ph_t_sq) const override;

	Real F_LL(constant::Hadron h, Real x, Real z, Real Q_sq, Real ph_t_sq) const override;
	Real F_LL_cos_phih(constant::Hadron h, Real x, Real z, Real Q_sq, Real ph_t_sq) const override;

	Real F_LT_cos_phih_m_phis(constant::Hadron h, Real x, Real z, Real Q_sq, Real ph_t_sq) const override;
	Real F_LT_cos_2phih_m_phis(constant::Hadron h, Real x, Real z, Real Q_sq, Real ph_t_sq) const override;
	Real F_LT_cos_phis(constant::Hadron h, Real x, Real z, Real Q_sq, Real ph_t_sq) const override;

	// Fragmentation functions.
	Real D1(constant::Hadron h, constant::Quark q, Real z, Real Q_sq) const;
	Real H1perpM1(constant::Hadron h, constant::Quark q, Real z, Real Q_sq) const;

	// Transverse momentum distributions.
	Real xf1(constant::Quark q, Real x, Real Q_sq) const;
	Real xf1TperpM1(constant::Quark q, Real x, Real Q_sq) const;
	Real xg1(constant::Quark q, Real x, Real Q_sq) const;
	Real xgT(constant::Quark q, Real x, Real Q_sq) const;
	Real xh1(constant::Quark q, Real x, Real Q_sq) const;
	Real xh1M1(constant::Quark q, Real x, Real Q_sq) const;
	Real xh1LperpM1(constant::Quark q, Real x, Real Q_sq) const;
	Real xh1TperpM2(constant::Quark q, Real x, Real Q_sq) const;
	Real xh1perpM1(constant::Quark q, Real x, Real Q_sq) const;
};

}
}
}

#endif

