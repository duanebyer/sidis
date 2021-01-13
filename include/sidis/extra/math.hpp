#ifndef SIDIS_MATH_HPP
#define SIDIS_MATH_HPP

namespace sidis {
namespace math {

inline float clamp(float min, float x, float max) {
	if (!(x > min)) {
		return min;
	} else if (!(x < max)) {
		return max;
	} else {
		return x;
	}
}

inline double clamp(double min, double x, double max) {
	if (!(x > min)) {
		return min;
	} else if (!(x < max)) {
		return max;
	} else {
		return x;
	}
}
inline long double clamp(long double min, long double x, long double max) {
	if (!(x > min)) {
		return min;
	} else if (!(x < max)) {
		return max;
	} else {
		return x;
	}
}

inline float sq(float x) {
	return x * x;
}
inline double sq(double x) {
	return x * x;
}
inline long double sq(long double x) {
	return x * x;
}

float sqrt1p_1m(float x);
double sqrt1p_1m(double x);
long double sqrt1p_1m(long double x);

float dilog(float x);
double dilog(double x);
long double dilog(long double x);

// These functions compute `N (k / K)` with reduced chance of overflow.
int prod_div(int N, int k, int K, int& rem);
long prod_div(long N, long k, long K, long& rem);
long long prod_div(long long N, long long k, long long K, long long& rem);

unsigned prod_div(unsigned N, unsigned k, unsigned K, unsigned& rem);
unsigned long prod_div(unsigned long N, unsigned long k, unsigned long K, unsigned long& rem);
unsigned long long prod_div(unsigned long long N, unsigned long long k, unsigned long long K, unsigned long long& rem);

inline int prod_div(int N, int k, int K) {
	int rem;
	return prod_div(N, k, K, rem);
}
inline long prod_div(long N, long k, long K) {
	long rem;
	return prod_div(N, k, K, rem);
}
inline long long prod_div(long long N, long long k, long long K) {
	long long rem;
	return prod_div(N, k, K, rem);
}

inline unsigned prod_div(unsigned N, unsigned k, unsigned K) {
	unsigned rem;
	return prod_div(N, k, K, rem);
}
inline unsigned long prod_div(unsigned long N, unsigned long k, unsigned long K) {
	unsigned long rem;
	return prod_div(N, k, K, rem);
}
inline unsigned long long prod_div(unsigned long long N, unsigned long long k, unsigned long long K) {
	unsigned long long rem;
	return prod_div(N, k, K, rem);
}

}
}

#endif

