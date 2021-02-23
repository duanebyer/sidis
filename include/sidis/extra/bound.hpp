#ifndef SIDIS_BOUNDS_HPP
#define SIDIS_BOUNDS_HPP

#include "sidis/numeric.hpp"

namespace sidis {
namespace math {

struct Bound {
	static Bound const INVALID;
	static Bound const ZERO;
	static Bound const UNIT;
	static Bound const FULL;
	static Bound const POSITIVE;
	static Bound const NEGATIVE;

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

	Bound() : _min(0.), _max(0.) { }
	Bound(Real min, Real max);
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
	bool operator==(Bound const& rhs) const {
		return _min == rhs._min && _max == rhs._max;
	}
	bool operator!=(Bound const& rhs) const {
		return !(*this == rhs);
	}

	// Combining bounds.
	Bound& operator&=(Bound const& rhs);
	Bound& operator|=(Bound const& rhs);

	// Arithmetic operations.
	Bound& operator+=(Real s);
	Bound& operator-=(Real s);
	Bound& operator*=(Real s);
	Bound& operator/=(Real s);
};

inline Bound operator&(Bound lhs, Bound const& rhs) {
	lhs &= rhs;
	return lhs;
}
inline Bound operator|(Bound lhs, Bound const& rhs) {
	lhs |= rhs;
	return lhs;
}

// Arithmetic operations.
inline Bound operator+(Bound lhs, Real rhs) {
	lhs += rhs;
	return lhs;
}
inline Bound operator+(Real lhs, Bound rhs) {
	rhs += lhs;
	return rhs;
}
inline Bound operator-(Bound lhs, Real rhs) {
	lhs -= rhs;
	return lhs;
}
inline Bound operator*(Bound lhs, Real rhs) {
	lhs *= rhs;
	return lhs;
}
inline Bound operator*(Real lhs, Bound rhs) {
	rhs *= lhs;
	return rhs;
}
inline Bound operator/(Bound lhs, Real rhs) {
	lhs /= rhs;
	return lhs;
}

}
}

#endif

