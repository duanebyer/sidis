#ifndef SIDIS_INTEGRATE_HPP
#define SIDIS_INTEGRATE_HPP

namespace sidis {
namespace math {

template<typename R, typename F>
R trapezoid(F f, R a, R b, unsigned n) {
	R delta = (b - a) / n;
	R result = (f(a) + f(b)) / 2;
	for (unsigned i = 1; i <= n - 1; ++i) {
		R x = i * (b - a) / n + a;
		result += f(x);
	}
	return result * delta;
}

}
}

#endif

