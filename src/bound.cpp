#include "sidis/bound.hpp"

#include <cmath>
#include <limits>

using namespace sidis;
using namespace sidis::math;

Bound::Bound() :
	_min(std::numeric_limits<Real>::quiet_NaN()),
	_max(std::numeric_limits<Real>::quiet_NaN()) {
}

Bound::Bound(Real min, Real max) : _min(min), _max(max) {
	if (!(_min <= _max)) {
		*this = BOUND_INVALID;
	}
}

Bound& Bound::operator&=(Bound const& rhs) {
	if (std::isnan(_min) || std::isnan(rhs._min)
			|| std::isnan(_max) || std::isnan(rhs._max)) {
		*this = BOUND_INVALID;
	} else {
		_min = std::fmax(_min, rhs._min);
		_max = std::fmin(_max, rhs._max);
		if (!(_min <= _max)) {
			*this = BOUND_INVALID;
		}
	}
	return *this;
}
Bound& Bound::operator|=(Bound const& rhs) {
	if (std::isnan(_min) || std::isnan(rhs._min)
			|| std::isnan(_max) || std::isnan(rhs._max)) {
		*this = BOUND_INVALID;
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
		*this = BOUND_INVALID;
	} else {
		_min *= s;
		_max *= s;
	}
	return *this;
}
Bound& Bound::operator/=(Real s) {
	if (s <= 0.) {
		*this = BOUND_INVALID;
	} else {
		_min /= s;
		_max /= s;
	}
	return *this;
}

