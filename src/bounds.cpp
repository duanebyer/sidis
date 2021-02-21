#include "sidis/extra/bounds.hpp"

#include <cmath>
#include <limits>

using namespace sidis;
using namespace sidis::math;

Bounds const Bounds::INVALID = Bounds(
	std::numeric_limits<Real>::quiet_NaN(),
	std::numeric_limits<Real>::quiet_NaN());
Bounds const Bounds::ZERO = Bounds(0., 0.);
Bounds const Bounds::UNIT = Bounds(0., 1.);
// TODO: Consider changing this to ±infinity.
Bounds const Bounds::FULL = Bounds(
	std::numeric_limits<Real>::lowest(),
	std::numeric_limits<Real>::max());
Bounds const Bounds::POSITIVE = Bounds(
	0.,
	std::numeric_limits<Real>::max());
Bounds const Bounds::NEGATIVE = Bounds(
	std::numeric_limits<Real>::lowest(),
	0.);

Bounds::Bounds(Real min, Real max) : _min(min), _max(max) {
	if (!(_min <= _max)) {
		*this = Bounds::INVALID;
	}
}

Bounds& Bounds::operator&=(Bounds const& rhs) {
	if (std::isnan(_min) || std::isnan(rhs._min)
			|| std::isnan(_max) || std::isnan(rhs._max)) {
		*this = Bounds::INVALID;
	} else {
		_min = std::fmax(_min, rhs._min);
		_max = std::fmin(_max, rhs._max);
		if (!(_min <= _max)) {
			*this = Bounds::INVALID;
		}
	}
	return *this;
}
Bounds& Bounds::operator|=(Bounds const& rhs) {
	if (std::isnan(_min) || std::isnan(rhs._min)
			|| std::isnan(_max) || std::isnan(rhs._max)) {
		*this = Bounds::INVALID;
	} else {
		_min = std::fmin(_min, rhs._min);
		_max = std::fmax(_max, rhs._max);
	}
	return *this;
}

Bounds& Bounds::operator+=(Real s) {
	_min += s;
	_max += s;
	return *this;
}
Bounds& Bounds::operator-=(Real s) {
	_min -= s;
	_max -= s;
	return *this;
}
Bounds& Bounds::operator*=(Real s) {
	if (s < 0.) {
		*this = Bounds::INVALID;
	} else {
		_min *= s;
		_max *= s;
	}
	return *this;
}
Bounds& Bounds::operator/=(Real s) {
	if (s <= 0.) {
		*this = Bounds::INVALID;
	} else {
		_min /= s;
		_max /= s;
	}
	return *this;
}

