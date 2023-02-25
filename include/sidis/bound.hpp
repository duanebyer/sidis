#ifndef SIDIS_BOUNDS_HPP
#define SIDIS_BOUNDS_HPP

#include "sidis/numeric.hpp"

namespace sidis {
namespace math {

/**
 * Half-open interval between two real numbers \f$[a, b)\f$.
 *
 * A Bound between \f$a\f$ and \f$b\f$ is required to satisfy \f$a \le b\f$. The
 * invalid Bound is also permitted through the use of BOUND_INVALID, and the
 * validity of the Bound can be checked through Bound::valid().
 */
class Bound {
private:
	Real _min;
	Real _max;

public:

	/// Lower end of the Bound.
	Real min() const {
		return _min;
	}
	/// Upper end of the Bound.
	Real max() const {
		return _max;
	}

	/// Initialize to BOUND_INVALID.
	Bound();
	/// Constructs a Bound between \p min and \p max. If the invariant
	/// `min <= max` is not satisfied, then sets the Bound to be the invalid
	/// Bound instead.
	Bound(Real min, Real max);
	/// Is \p x contained in the interval?
	bool contains(Real x) const {
		return _min <= x && x <= _max;
	}
	/// Are both ends of \p bound contained within this Bound?
	bool contains(Bound other) const {
		return contains(other._min) && contains(other._max);
	}
	/// Linearly interpolate between the ends of the Bound with \p x in
	/// \f$[0, 1)\f$.
	Real lerp(Real x) const {
		return (1. - x) * _min + x * _max;
	}
	/// Return whether the Bound is valid.
	bool valid() const {
		return _min <= _max;
	}
	/// The distance between the ends of the Bound.
	Real size() const {
		return _max - _min;
	}

	/// Checks for equality. Note that an invalid Bound is not equal to itself.
	bool operator==(Bound const& rhs) const {
		return _min == rhs._min && _max == rhs._max;
	}
	/// Checks for inequality. Note that an invalid Bound is not equal to
	/// itself.
	bool operator!=(Bound const& rhs) const {
		return !(*this == rhs);
	}

	/// \copydoc Bound::operator&(Bound, Bound const&)
	Bound& operator&=(Bound const& rhs);
	/// \copydoc Bound::operator|(Bound, Bound const&)
	Bound& operator|=(Bound const& rhs);

	/// \copydoc Bound::operator+(Bound, Real)
	Bound& operator+=(Real s);
	/// \copydoc Bound::operator-(Bound, Real)
	Bound& operator-=(Real s);
	/// \copydoc Bound::operator*(Bound, Real)
	Bound& operator*=(Real s);
	/// \copydoc Bound::operator/(Bound, Real)
	Bound& operator/=(Real s);
};

/// Unique invalid Bound.
Bound const BOUND_INVALID;
/// Bound for \f$[0, 0)\f$.
Bound const BOUND_ZERO(0., 0.);
/// Bound for \f$[0, 1)\f$.
Bound const BOUND_UNIT(0., 1.);
/// Bound for \f$[-\infty, +\infty)\f$.
Bound const BOUND_FULL(-1. / 0., 1. / 0.);
/// Bound for \f$[0, +\infty)\f$.
Bound const BOUND_POSITIVE(0., 1. / 0.);
/// Bound for \f$[-\infty, 0)\f$.
Bound const BOUND_NEGATIVE(-1. / 0., 0.);

/// Intersection of two Bound%s. If there is no overlap, then give
/// BOUND_INVALID.
/// \relates Bound
inline Bound operator&(Bound lhs, Bound const& rhs) {
	lhs &= rhs;
	return lhs;
}
/// Union of two Bound%s. Includes any region in-between the two Bound%s that
/// belongs to neither of them.
/// \relates Bound
inline Bound operator|(Bound lhs, Bound const& rhs) {
	lhs |= rhs;
	return lhs;
}

/// Shifts the endpoints of the Bound.
/// \relates Bound
inline Bound operator+(Bound lhs, Real rhs) {
	lhs += rhs;
	return lhs;
}
/// \copydoc Bound::operator+(Bound, Real)
/// \relates Bound
inline Bound operator+(Real lhs, Bound rhs) {
	rhs += lhs;
	return rhs;
}
/// \copydoc Bound::operator+(Bound, Real)
/// \relates Bound
inline Bound operator-(Bound lhs, Real rhs) {
	lhs -= rhs;
	return lhs;
}

/// Scales both endpoints of the Bound. Scaling by a negative number results in
/// an invalid Bound.
/// \relates Bound
inline Bound operator*(Bound lhs, Real rhs) {
	lhs *= rhs;
	return lhs;
}
/// \copydoc Bound::operator*(Bound, Real)
/// \relates Bound
inline Bound operator*(Real lhs, Bound rhs) {
	rhs *= lhs;
	return rhs;
}
/// \copydoc Bound::operator*(Bound, Real)
/// \relates Bound
inline Bound operator/(Bound lhs, Real rhs) {
	lhs /= rhs;
	return lhs;
}

}
}

#endif

