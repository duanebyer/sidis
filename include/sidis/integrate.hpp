#ifndef SIDIS_INTEGRATE_HPP
#define SIDIS_INTEGRATE_HPP

namespace sidis {
namespace math {

template<typename R, typename F>
R reimann(F f, R a, R b, unsigned n) {
	R delta = (b - a) / n;
	R result = 0.;
	for (unsigned i = 0; i < n; ++i) {
		R x = (i * b + (n - i) * a + 0.5) / n;
		result += f(x);
	}
	return result * delta;
}

template<typename R, typename F>
R trapezoid(F f, R a, R b, unsigned n) {
	R delta = (b - a) / n;
	R result = 0.5 * (f(a) + f(b));
	for (unsigned i = 1; i <= n - 1; ++i) {
		R x = (i * b + (n - i) * a) / n;
		result += f(x);
	}
	return result * delta;
}

}
}

#endif

