#include "sidis/math.hpp"

#include <cmath>
#include <limits>

using namespace sidis;
using namespace sidis::math;

template<typename T>
static T dilog_impl(T x) {
	// Ensure we have as many digits of pi as can fit the `T` type.
	T pi = 3.141592653589793238462643383279502797479068098137295573004504331874296718662975536062731407582759857177734375;

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
		b = -sq(pi) / 6. - 0.5 * sq(std::log(x));
		x = 1. / x;
	} else if (x > 1.) {
		a = 1.;
		b = sq(pi) / 6. + 0.5 * std::log(x) * std::log(x / sq(1. - x));
		x = (x - 1.) / x;
	} else if (x > 0.5) {
		a = -1.;
		b = sq(pi) / 6. - std::log(x) * std::log(1. - x);
		x = 1. - x;
	} else if (x > -0.5) {
		a = 1.;
		b = 0.;
		x = x;
	} else if (x > -1.) {
		a = -1.;
		b = -0.5 * sq(std::log(1. - x));
		x = x / (x - 1.);
	} else {
		a = 1.;
		b = -sq(pi) / 6. + 0.5 * std::log(1. - x) * std::log((1. - x) / sq(x));
		x = 1. / (1. - x);
	}

	// Size of the mantissa in base 2.
	double d = std::numeric_limits<T>::digits
		* std::log2(std::numeric_limits<T>::radix);
	// Number of terms needed to reach desired precision.
	unsigned n_max = (unsigned) std::ceil(
		(0.4 * d + 6. * (1. - std::log(d))) * x + 0.2 * d);
	T result = 0.;
	T numerator = x;
	for (unsigned n = 1; n < n_max; ++n) {
		result += numerator / (n * n);
		numerator *= x;
	}

	return a * result + b;
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

