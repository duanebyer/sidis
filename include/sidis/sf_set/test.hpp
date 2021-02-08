#ifndef SIDIS_SF_SET_TEST_HPP
#define SIDIS_SF_SET_TEST_HPP

#include "sidis/structure_function.hpp"

namespace sidis {
namespace sf {
namespace set {

/*
 * Structure function set that returns a constant value of one for some subset
 * of all structure functions, for testing purposes.
 *
 * This class takes a bitset determining which of the structure functions will
 * return a constant value of one, and which will return zero. This is used for
 * doing tests on the contributions of structure functions which normally might
 * not significantly affect the cross-section.
 */
class TestSfSet final : public SfSet {
	bool mask[18];

public:
	TestSfSet(constant::Nucleus target, const bool (&mask)[18]);

	Real F_UUL(constant::Hadron h, Real x, Real z, Real Q_sq, Real ph_t_sq) const override;
	Real F_UUT(constant::Hadron h, Real x, Real z, Real Q_sq, Real ph_t_sq) const override;
	Real F_UU_cos_phih(constant::Hadron h, Real x, Real z, Real Q_sq, Real ph_t_sq) const override;
	Real F_UU_cos_2phih(constant::Hadron h, Real x, Real z, Real Q_sq, Real ph_t_sq) const override;

	Real F_UL_sin_phih(constant::Hadron h, Real x, Real z, Real Q_sq, Real ph_t_sq) const override;
	Real F_UL_sin_2phih(constant::Hadron h, Real x, Real z, Real Q_sq, Real ph_t_sq) const override;

	Real F_UTL_sin_phih_m_phis(constant::Hadron h, Real x, Real z, Real Q_sq, Real ph_t_sq) const override;
	Real F_UTT_sin_phih_m_phis(constant::Hadron h, Real x, Real z, Real Q_sq, Real ph_t_sq) const override;
	Real F_UT_sin_2phih_m_phis(constant::Hadron h, Real x, Real z, Real Q_sq, Real ph_t_sq) const override;
	Real F_UT_sin_3phih_m_phis(constant::Hadron h, Real x, Real z, Real Q_sq, Real ph_t_sq) const override;
	Real F_UT_sin_phis(constant::Hadron h, Real x, Real z, Real Q_sq, Real ph_t_sq) const;
	Real F_UT_sin_phih_p_phis(constant::Hadron h, Real x, Real z, Real Q_sq, Real ph_t_sq) const override;

	Real F_LU_sin_phih(constant::Hadron h, Real x, Real z, Real Q_sq, Real ph_t_sq) const override;

	Real F_LL(constant::Hadron h, Real x, Real z, Real Q_sq, Real ph_t_sq) const override;
	Real F_LL_cos_phih(constant::Hadron h, Real x, Real z, Real Q_sq, Real ph_t_sq) const override;

	Real F_LT_cos_phih_m_phis(constant::Hadron h, Real x, Real z, Real Q_sq, Real ph_t_sq) const override;
	Real F_LT_cos_2phih_m_phis(constant::Hadron h, Real x, Real z, Real Q_sq, Real ph_t_sq) const override;
	Real F_LT_cos_phis(constant::Hadron h, Real x, Real z, Real Q_sq, Real ph_t_sq) const override;
};

}
}
}

#endif

