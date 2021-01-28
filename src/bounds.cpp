#include "sidis/extra/bounds.hpp"

#include <cmath>
#include <limits>

using namespace sidis;
using namespace sidis::math;

Bounds const Bounds::INVALID = Bounds(
	std::numeric_limits<Real>::quiet_NaN(),
	std::numeric_limits<Real>::quiet_NaN());
Bounds const Bounds::FULL = Bounds(
	std::numeric_limits<Real>::lowest(),
	std::numeric_limits<Real>::max());

Bounds math::operator&(Bounds lhs, Bounds rhs) {
	return Bounds(
		std::fmax(lhs.min, rhs.min),
		std::fmin(lhs.max, rhs.max));
}

Bounds math::operator|(Bounds lhs, Bounds rhs) {
	return Bounds(
		std::fmin(lhs.min, rhs.min),
		std::fmax(lhs.max, rhs.max));
}

