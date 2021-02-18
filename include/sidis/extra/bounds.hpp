#ifndef SIDIS_BOUNDS_HPP
#define SIDIS_BOUNDS_HPP

#include "sidis/numeric.hpp"

namespace sidis {
namespace math {

struct Bounds {
	static Bounds const INVALID;
	static Bounds const ZERO;
	static Bounds const UNIT;
	static Bounds const FULL;
	static Bounds const POSITIVE;
	static Bounds const NEGATIVE;

private:
	Real _min;
	Real _max;

public:

	Real min() const {
		return _min;
	}
	Real max() const {
		return _max;
	}

	Bounds(Real min, Real max);
	bool contains(Real x) const {
		return _min <= x && x < _max;
	}
	Real lerp(Real x) const {
		return (1. - x) * _min + x * _max;
	}
	Real valid() const {
		return _min <= _max;
	}
	Real size() const {
		return _max - _min;
	}

	// Equality operations.
	bool operator==(Bounds const& rhs) const {
		return _min == rhs._min && _max == rhs._max;
	}
	bool operator!=(Bounds const& rhs) const {
		return !(*this == rhs);
	}

	// Combining bounds.
	Bounds& operator&=(Bounds const& rhs);
	Bounds& operator|=(Bounds const& rhs);

	// Arithmetic operations.
	Bounds& operator+=(Real s);
	Bounds& operator-=(Real s);
	Bounds& operator*=(Real s);
	Bounds& operator/=(Real s);
};

inline Bounds operator&(Bounds lhs, Bounds const& rhs) {
	lhs &= rhs;
	return lhs;
}
inline Bounds operator|(Bounds lhs, Bounds const& rhs) {
	lhs |= rhs;
	return lhs;
}

// Arithmetic operations.
inline Bounds operator+(Bounds lhs, Real rhs) {
	lhs += rhs;
	return lhs;
}
inline Bounds operator+(Real lhs, Bounds rhs) {
	rhs += lhs;
	return rhs;
}
inline Bounds operator-(Bounds lhs, Real rhs) {
	lhs -= rhs;
	return lhs;
}
inline Bounds operator*(Bounds lhs, Real rhs) {
	lhs *= rhs;
	return lhs;
}
inline Bounds operator*(Real lhs, Bounds rhs) {
	rhs *= lhs;
	return rhs;
}
inline Bounds operator/(Bounds lhs, Real rhs) {
	lhs /= rhs;
	return lhs;
}

}
}

#endif

