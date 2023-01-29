#include "sidis/extra/math.hpp"

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
	T c = 1.64493406684822643647241516664602518921894990120679843773555822937000\
74704032008738336289006197587053040043189623371906796287247L;

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
		//x = x;
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

template<typename T>
T prod_div_impl(T N, T k, T K, T& rem) {
	if (k < 0) {
		N = -N;
		k = -k;
	}
	T N_max = N >= 0 ?
		std::numeric_limits<T>::max() :
		std::numeric_limits<T>::min();
	if (k < N_max / N) {
		// Base case, when there is no risk of overflow.
		rem = (k * N) % K;
		return (k * N) / K;
	} else {
		// The basic algorithm used here is to take `N / K`, and then scale
		// repeatedly by fractions of the form `(p + 1) / p`, for `p = 1..k`.
		T k_0 = N_max / N;
		T quot = (N * k_0) / K;
		rem = (N * k_0) % K;
		for (T k_p = k_0; k_p < k; ++k_p) {
			T quot_p = quot / k_p;
			T rem_p = quot % k_p;
			T prod_num = (rem * (k_p + 1) + rem_p * K);
			T prod = prod_num / k_p;
			quot += quot_p + prod / K;
			rem = prod % K;
		}
		return quot;
	}
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

int math::prod_div(int N, int k, int K, int& rem) {
	return prod_div_impl<int>(N, k, K, rem);
}
long math::prod_div(long N, long k, long K, long& rem) {
	return prod_div_impl<long>(N, k, K, rem);
}
long long math::prod_div(long long N, long long k, long long K, long long& rem) {
	return prod_div_impl<long long>(N, k, K, rem);
}

unsigned math::prod_div(unsigned N, unsigned k, unsigned K, unsigned& rem) {
	return prod_div_impl<unsigned>(N, k, K, rem);
}
unsigned long math::prod_div(unsigned long N, unsigned long k, unsigned long K, unsigned long& rem) {
	return prod_div_impl<unsigned long>(N, k, K, rem);
}
unsigned long long math::prod_div(unsigned long long N, unsigned long long k, unsigned long long K, unsigned long long& rem) {
	return prod_div_impl<unsigned long long>(N, k, K, rem);
}

