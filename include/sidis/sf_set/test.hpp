#ifndef SIDIS_SF_SET_TEST_HPP
#define SIDIS_SF_SET_TEST_HPP

#include "sidis/structure_function.hpp"

namespace sidis {
namespace sf {
namespace set {

/**
 * Structure function set that always returns 1, for testing purposes.
 */
class TestSfSet final : public SfSet {
public:
	/// Constructs an SfSet for \p target that returns 1 for every structure
	/// function.
	TestSfSet(part::Nucleus target) : SfSet(target) { }
	TestSfSet(TestSfSet const& other) : TestSfSet(other.target) { }
	TestSfSet(TestSfSet&& other) : TestSfSet(other.target) { }
	TestSfSet& operator=(TestSfSet const& other) = delete;
	TestSfSet& operator=(TestSfSet&& other) = delete;
	virtual ~TestSfSet() = default;

	Real F_UUL(part::Hadron, Real, Real, Real, Real) const override {
		return 1.;
	}
	Real F_UUT(part::Hadron, Real, Real, Real, Real) const override {
		return 1.;
	}
	Real F_UU_cos_phih(part::Hadron, Real, Real, Real, Real) const override {
		return 1.;
	}
	Real F_UU_cos_2phih(part::Hadron, Real, Real, Real, Real) const override {
		return 1.;
	}

	Real F_UL_sin_phih(part::Hadron, Real, Real, Real, Real) const override {
		return 1.;
	}
	Real F_UL_sin_2phih(part::Hadron, Real, Real, Real, Real) const override {
		return 1.;
	}

	Real F_UTL_sin_phih_m_phis(part::Hadron, Real, Real, Real, Real) const override {
		return 1.;
	}
	Real F_UTT_sin_phih_m_phis(part::Hadron, Real, Real, Real, Real) const override {
		return 1.;
	}
	Real F_UT_sin_2phih_m_phis(part::Hadron, Real, Real, Real, Real) const override {
		return 1.;
	}
	Real F_UT_sin_3phih_m_phis(part::Hadron, Real, Real, Real, Real) const override {
		return 1.;
	}
	Real F_UT_sin_phis(part::Hadron, Real, Real, Real, Real) const override {
		return 1.;
	}
	Real F_UT_sin_phih_p_phis(part::Hadron, Real, Real, Real, Real) const override {
		return 1.;
	}

	Real F_LU_sin_phih(part::Hadron, Real, Real, Real, Real) const override {
		return 1.;
	}

	Real F_LL(part::Hadron, Real, Real, Real, Real) const override {
		return 1.;
	}
	Real F_LL_cos_phih(part::Hadron, Real, Real, Real, Real) const override {
		return 1.;
	}

	Real F_LT_cos_phih_m_phis(part::Hadron, Real, Real, Real, Real) const override {
		return 1.;
	}
	Real F_LT_cos_2phih_m_phis(part::Hadron, Real, Real, Real, Real) const override {
		return 1.;
	}
	Real F_LT_cos_phis(part::Hadron, Real, Real, Real, Real) const override {
		return 1.;
	}
};

}
}
}

#endif

