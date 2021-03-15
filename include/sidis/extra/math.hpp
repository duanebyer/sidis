#ifndef SIDIS_MATH_HPP
#define SIDIS_MATH_HPP

namespace sidis {
namespace math {

/**
 * \defgroup MathGroup Math functions
 * Useful math functions.
 */
/// \{

/// Clamp a number between an upper and lower bound.
inline float clamp(float min, float x, float max) {
	if (!(x > min)) {
		return min;
	} else if (!(x < max)) {
		return max;
	} else {
		return x;
	}
}
/// \copydoc clamp()
inline double clamp(double min, double x, double max) {
	if (!(x > min)) {
		return min;
	} else if (!(x < max)) {
		return max;
	} else {
		return x;
	}
}
/// \copydoc clamp()
inline long double clamp(long double min, long double x, long double max) {
	if (!(x > min)) {
		return min;
	} else if (!(x < max)) {
		return max;
	} else {
		return x;
	}
}

/// Square a number.
inline float sq(float x) {
	return x * x;
}
/// \copydoc sq()
inline double sq(double x) {
	return x * x;
}
/// \copydoc sq()
inline long double sq(long double x) {
	return x * x;
}

/// Compute \f$\sqrt{1 + x} - 1\f$ with better numerical precision.
float sqrt1p_1m(float x);
/// \copydoc sqrt1p_1m()
double sqrt1p_1m(double x);
/// \copydoc sqrt1p_1m()
long double sqrt1p_1m(long double x);

/// Compute the dilog function with good numerical precision.
float dilog(float x);
/// \copydoc dilog()
double dilog(double x);
/// \copydoc dilog()
long double dilog(long double x);

/// Compute \f$N (k / K)\f$ for integers with reduced chance of overflow.
int prod_div(int N, int k, int K, int& rem);
/// \copydoc prod_div()
long prod_div(long N, long k, long K, long& rem);
/// \copydoc prod_div()
long long prod_div(long long N, long long k, long long K, long long& rem);

/// \copydoc prod_div()
unsigned prod_div(unsigned N, unsigned k, unsigned K, unsigned& rem);
/// \copydoc prod_div()
unsigned long prod_div(unsigned long N, unsigned long k, unsigned long K, unsigned long& rem);
/// \copydoc prod_div()
unsigned long long prod_div(unsigned long long N, unsigned long long k, unsigned long long K, unsigned long long& rem);

/// \copydoc prod_div()
inline int prod_div(int N, int k, int K) {
	int rem;
	return prod_div(N, k, K, rem);
}
/// \copydoc prod_div()
inline long prod_div(long N, long k, long K) {
	long rem;
	return prod_div(N, k, K, rem);
}
/// \copydoc prod_div()
inline long long prod_div(long long N, long long k, long long K) {
	long long rem;
	return prod_div(N, k, K, rem);
}

/// \copydoc prod_div()
inline unsigned prod_div(unsigned N, unsigned k, unsigned K) {
	unsigned rem;
	return prod_div(N, k, K, rem);
}
/// \copydoc prod_div()
inline unsigned long prod_div(unsigned long N, unsigned long k, unsigned long K) {
	unsigned long rem;
	return prod_div(N, k, K, rem);
}
/// \copydoc prod_div()
inline unsigned long long prod_div(unsigned long long N, unsigned long long k, unsigned long long K) {
	unsigned long long rem;
	return prod_div(N, k, K, rem);
}
/// \}

}
}

#endif

