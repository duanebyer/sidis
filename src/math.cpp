#include "sidis/math.hpp"

#include <cmath>
#include <limits>

using namespace sidis;
using namespace sidis::math;

namespace {

template<typename T>
T sqrt1p_1m_impl(T x) {
	return std::expm1(0.5*std::log1p(x));
}

// Compute `log(|x|) log(|1 - x|)`
template<typename T>
T log_log1m(T x) {
	if (!std::isfinite(x)) {
		if (std::isinf(x)) {
			return std::numeric_limits<T>::infinity();
		} else {
			return x;
		}
	} else if (x > 1.) {
		return std::log(x) * std::log(x - 1.);
	} else if (x == 1.) {
		return 0.;
	} else if (x > 0.5) {
		return std::log1p(x - 1.) * std::log(1. - x);
	} else if (x > 0.) {
		return std::log(x) * std::log1p(-x);
	} else {
		return std::log(-x) * std::log1p(-x);
	}
}

template<typename T>
T dilog_impl(T x) {
	// Ensure we have as many digits of `π^2/6` as can fit the `T` type, up to
	// `long double`.
	T c = 1.644934066848226436472415166646025189219L;

	// To compute the dilogarithm, first bring the argument into the range
	// (-0.5, 0.5).
	T a = 1.;
	T b = 0.;
	if (!std::isfinite(x)) {
		if (std::isinf(x)) {
			return -std::numeric_limits<T>::infinity();
		} else {
			return x;
		}
	} else if (x > 2.) {
		a = -1.;
		b = 2. * c - 0.5 * sq(std::log(x));
		x = 1. / x;
	} else if (x > 1.) {
		a = 1.;
		b = c + 0.5 * sq(std::log(x)) - log_log1m(x);
		x = (x - 1.) / x;
	} else if (x > 0.5) {
		a = -1.;
		b = c - log_log1m(x);
		x = 1. - x;
	} else if (x > -0.5) {
		a = 1.;
		b = 0.;
		x = x;
	} else if (x > -1.) {
		a = -1.;
		b = -0.5 * sq(std::log1p(-x));
		x = x / (x - 1.);
	} else {
		a = 1.;
		b = -c + 0.5 * sq(std::log1p(-x)) - log_log1m(x);
		x = 1. / (1. - x);
	}

	// Size of the mantissa in base 2.
	double d = std::numeric_limits<T>::digits
		* std::log2(std::numeric_limits<T>::radix);
	// Number of terms needed to reach desired precision. This bound was
	// experimentally determined to satisfy this condition for `d` greater than
	// 10 and `x` in range (-0.5, 0.5).
	unsigned n_max = (unsigned) std::ceil(
		(1.4 * d + 6. * (1. - std::log(d))) * std::abs(x) + 0.3 * d) + 1;
	T result = 0.;
	T numerator = x;
	for (unsigned n = 1; n < n_max; ++n) {
		result += numerator / (n * n);
		numerator *= x;
	}

	return a * result + b;
}

}

float math::sqrt1p_1m(float x) {
	return sqrt1p_1m_impl<float>(x);
}
double math::sqrt1p_1m(double x) {
	return sqrt1p_1m_impl<double>(x);
}
long double math::sqrt1p_1m(long double x) {
	return sqrt1p_1m_impl<long double>(x);
}

float math::dilog(float x) {
	return dilog_impl<float>(x);
}
double math::dilog(double x) {
	return dilog_impl<double>(x);
}
long double math::dilog(long double x) {
	return dilog_impl<long double>(x);
}

