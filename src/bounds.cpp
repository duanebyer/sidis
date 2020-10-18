#include "sidis/extra/bounds.hpp"

#include <cmath>

using namespace sidis;
using namespace sidis::math;

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

