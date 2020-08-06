#ifndef SIDIS_MATH_HPP
#define SIDIS_MATH_HPP

namespace sidis {
namespace math {

inline float sq(float x) {
	return x * x;
}
inline double sq(double x) {
	return x * x;
}
inline long double sq(long double x) {
	return x * x;
}

inline float lerp(float lower, float upper, float x) {
	return (1. - x) * lower + x * upper;
}
inline double lerp(double lower, double upper, double x) {
	return (1. - x) * lower + x * upper;
}
inline long double lerp(long double lower, long double upper, long double x) {
	return (1. - x) * lower + x * upper;
}

float dilog(float x);
double dilog(double x);
long double dilog(long double x);

}
}

#endif

