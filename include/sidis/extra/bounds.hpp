#ifndef SIDIS_BOUNDS_HPP
#define SIDIS_BOUNDS_HPP

#include "sidis/numeric.hpp"

namespace sidis {
namespace math {

struct Bounds {
	static Bounds const INVALID;
	static Bounds const FULL;

	Real min;
	Real max;
	Bounds(Real min, Real max) : min(min), max(max) { }
	bool operator()(Real x) const {
		return min <= x && x < max;
	}
	Real valid() const {
		return min < max;
	}
	Real lerp(Real x) const {
		return (1. - x) * min + x * max;
	}
	Real size() const {
		return max - min;
	}
};

Bounds operator&(Bounds lhs, Bounds rhs);
Bounds operator|(Bounds lhs, Bounds rhs);

}
}

#endif

