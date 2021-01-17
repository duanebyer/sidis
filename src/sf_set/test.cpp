#include "sidis/sf_set/test.hpp"

#include <cstddef>

using namespace sidis;
using namespace sidis::constant;
using namespace sidis::sf;
using namespace sidis::sf::model;

TestSfSet::TestSfSet(Nucleus target, const bool (&mask)[18]) :
		SfSet(target),
		mask() {
	for (std::size_t idx = 0; idx < sizeof(mask) / sizeof(bool); ++idx) {
		this->mask[idx] = mask[idx];
	}
}

Real TestSfSet::F_UUL(constant::Hadron, Real, Real, Real, Real) const {
	return static_cast<Real>(mask[0]);
}
Real TestSfSet::F_UUT(constant::Hadron, Real, Real, Real, Real) const {
	return static_cast<Real>(mask[1]);
}
Real TestSfSet::F_UU_cos_phih(constant::Hadron, Real, Real, Real, Real) const {
	return static_cast<Real>(mask[2]);
}
Real TestSfSet::F_UU_cos_2phih(constant::Hadron, Real, Real, Real, Real) const {
	return static_cast<Real>(mask[3]);
}

Real TestSfSet::F_UL_sin_phih(constant::Hadron, Real, Real, Real, Real) const {
	return static_cast<Real>(mask[4]);
}
Real TestSfSet::F_UL_sin_2phih(constant::Hadron, Real, Real, Real, Real) const {
	return static_cast<Real>(mask[5]);
}

Real TestSfSet::F_UTL_sin_phih_m_phis(constant::Hadron, Real, Real, Real, Real) const {
	return static_cast<Real>(mask[6]);
}
Real TestSfSet::F_UTT_sin_phih_m_phis(constant::Hadron, Real, Real, Real, Real) const {
	return static_cast<Real>(mask[7]);
}
Real TestSfSet::F_UT_sin_2phih_m_phis(constant::Hadron, Real, Real, Real, Real) const {
	return static_cast<Real>(mask[8]);
}
Real TestSfSet::F_UT_sin_3phih_m_phis(constant::Hadron, Real, Real, Real, Real) const {
	return static_cast<Real>(mask[9]);
}
Real TestSfSet::F_UT_sin_phis(constant::Hadron, Real, Real, Real, Real) const {
	return static_cast<Real>(mask[10]);
}
Real TestSfSet::F_UT_sin_phih_p_phis(constant::Hadron, Real, Real, Real, Real) const {
	return static_cast<Real>(mask[11]);
}

Real TestSfSet::F_LU_sin_phih(constant::Hadron, Real, Real, Real, Real) const {
	return static_cast<Real>(mask[12]);
}

Real TestSfSet::F_LL(constant::Hadron, Real, Real, Real, Real) const {
	return static_cast<Real>(mask[13]);
}
Real TestSfSet::F_LL_cos_phih(constant::Hadron, Real, Real, Real, Real) const {
	return static_cast<Real>(mask[14]);
}

Real TestSfSet::F_LT_cos_phih_m_phis(constant::Hadron, Real, Real, Real, Real) const {
	return static_cast<Real>(mask[15]);
}
Real TestSfSet::F_LT_cos_2phih_m_phis(constant::Hadron, Real, Real, Real, Real) const {
	return static_cast<Real>(mask[16]);
}
Real TestSfSet::F_LT_cos_phis(constant::Hadron, Real, Real, Real, Real) const {
	return static_cast<Real>(mask[17]);
}

