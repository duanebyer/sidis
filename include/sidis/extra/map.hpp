#ifndef SIDIS_MAP_HPP
#define SIDIS_MAP_HPP

#include <sidis/numeric.hpp>
#include <sidis/bound.hpp>

namespace sidis {
namespace math {
namespace map {

// Trivial linear transformation.
struct Linear {
	Real u(Real x) const {
		return x;
	}
	Real x(Real u, Real* jac) const {
		*jac = 1.;
		return u;
	}
};

// Transforms according to reciprocal function.
struct Inverse {
	Real offset;
	Inverse(Real offset=0.) : offset(offset) { }
	Real u(Real x) const {
		return -1. / (x + offset);
	}
	Real x(Real u, Real* jac) const {
		*jac = sq(1. / u);
		return -1. / u - offset;
	}
};

// Transforms according to logarithm.
struct Log {
	Real offset;
	Log(Real offset=0.) : offset(offset) { }
	Real u(Real x) const {
		return std::log(x + offset);
	}
	Real x(Real u, Real* jac) const {
		Real exp = std::exp(u);
		*jac = exp;
		return exp - offset;
	}
};

// Transforms according to decaying exponential.
struct Decay {
	Real length;
	Decay(Real length=1.) : length(length) { }
	Real u(Real x) const {
		return -std::expm1(-x / length);
	}
	Real x(Real u, Real* jac) const {
		*jac = length / (1. - u);
		return -length * std::log1p(-u);
	}
};

// Transforms according to sigmoid. Useful for peaked functions.
struct Sigmoid {
	Real center;
	Real width;
	Sigmoid(Real center=0., Real width=1.) : center(center), width(width) { }
	Real u(Real x) const {
		// Use the `arcsinh` function here, as it has logarithmic tails which
		// allows for the non-peak part of the integrand to still be captured
		// accurately, even with fat tails.
		return std::asinh((x - center) / width);
	}
	Real x(Real u, Real* jac) const {
		*jac = width * std::cosh(u);
		return width * std::sinh(u) + center;
	}
};

// Transforms according to a double-sigmoid. Useful for doubly-peaked
// integrands. One of the peaks is specified in the regular way, and then the
// other peak is specified by the ratio of its center and width to those of the
// first peak.
struct Sigmoid2 {
	Real center;
	Real width;
	// Ratio between second sigmoid center and first one.
	Real center_r;
	// Ratio between second sigmoid width and first one.
	Real width_r;
	Sigmoid2(Real center=0., Real width=1., Real center_r=-1., Real width_r=1.) :
		center(center),
		width(width),
		center_r(center_r),
		width_r(width_r) { }
	Real u(Real x) const {
		Real t1 = (x - center) / width;
		Real t2 = (x - center * center_r) / (width * width_r);
		return std::asinh(t1) + std::asinh(t2);
	}
	Real x(Real u, Real* jac) const {
		// Against all odds, this function really can be inverted.
		Real sinhu = std::sinh(u);
		Real coshu = std::cosh(u);
		Real a = center_r + sq(width_r) + width_r * (1. + center_r) * coshu;
		Real b = 1. + sq(width_r) + 2. * width_r * coshu;
		Real c = sq((center_r - 1.) * (center / width));
		Real det = std::sqrt(b + c);
		Real a_p = width_r * (1. + center_r) * sinhu;
		Real b_p = 2. * width_r * sinhu;
		Real num = center * a + width * width_r * det * sinhu;
		*jac = (center * a_p
			+ width * width_r * (det * coshu + 0.5 * sinhu * b_p / det)
			- num * b_p / b) / b;
		return num / b;
	}
};

}

// Uses the specified transform to take a point `p` in the range `[0, 1)` and
// map it into another `bound`, with a resulting Jacobian `jac`.
template<typename T>
Real apply_map(T const& t, Real p, Bound bound, Real* jac) {
	Bound u_bound(t.u(bound.min()), t.u(bound.max()));
	Real u = u_bound.lerp(p);
	Real x = t.x(u, jac);
	*jac *= u_bound.size();
	return x;
}

}
}

#endif

