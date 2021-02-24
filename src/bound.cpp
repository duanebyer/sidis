#include "sidis/bound.hpp"

#include <cmath>
#include <limits>

using namespace sidis;
using namespace sidis::math;

Bound const Bound::INVALID = Bound(
	std::numeric_limits<Real>::quiet_NaN(),
	std::numeric_limits<Real>::quiet_NaN());
Bound const Bound::ZERO = Bound(0., 0.);
Bound const Bound::UNIT = Bound(0., 1.);
Bound const Bound::FULL = Bound(
	-std::numeric_limits<Real>::infinity(),
	std::numeric_limits<Real>::infinity());
Bound const Bound::POSITIVE = Bound(
	0.,
	std::numeric_limits<Real>::max());
Bound const Bound::NEGATIVE = Bound(
	std::numeric_limits<Real>::lowest(),
	0.);

Bound::Bound(Real min, Real max) : _min(min), _max(max) {
	if (!(_min <= _max)) {
		*this = Bound::INVALID;
	}
}

Bound& Bound::operator&=(Bound const& rhs) {
	if (std::isnan(_min) || std::isnan(rhs._min)
			|| std::isnan(_max) || std::isnan(rhs._max)) {
		*this = Bound::INVALID;
	} else {
		_min = std::fmax(_min, rhs._min);
		_max = std::fmin(_max, rhs._max);
		if (!(_min <= _max)) {
			*this = Bound::INVALID;
		}
	}
	return *this;
}
Bound& Bound::operator|=(Bound const& rhs) {
	if (std::isnan(_min) || std::isnan(rhs._min)
			|| std::isnan(_max) || std::isnan(rhs._max)) {
		*this = Bound::INVALID;
	} else {
		_min = std::fmin(_min, rhs._min);
		_max = std::fmax(_max, rhs._max);
	}
	return *this;
}

Bound& Bound::operator+=(Real s) {
	_min += s;
	_max += s;
	return *this;
}
Bound& Bound::operator-=(Real s) {
	_min -= s;
	_max -= s;
	return *this;
}
Bound& Bound::operator*=(Real s) {
	if (s < 0.) {
		*this = Bound::INVALID;
	} else {
		_min *= s;
		_max *= s;
	}
	return *this;
}
Bound& Bound::operator/=(Real s) {
	if (s <= 0.) {
		*this = Bound::INVALID;
	} else {
		_min /= s;
		_max /= s;
	}
	return *this;
}

