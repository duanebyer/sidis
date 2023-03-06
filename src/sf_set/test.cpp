#include "sidis/sf_set/test.hpp"

#include <cmath>

#include "sidis/extra/math.hpp"

using namespace sidis;
using namespace sidis::math;
using namespace sidis::sf;
using namespace sidis::sf::set;

namespace {

// TODO: Put this in the math header, and implement it more accurately as alpha
// and beta go to zero.
Real L(Real x, Real alpha, Real beta) {
	if (alpha == 0.) {
		return std::pow(1. - x, beta);
	} else if (beta == 0.) {
		return std::pow(x, alpha);
	} else {
		Real norm = std::pow(std::abs(alpha) + std::abs(beta), alpha + beta)
			/ (std::pow(std::abs(alpha), alpha) * std::pow(std::abs(beta), beta));
		Real scale = std::pow(x, alpha) * std::pow(1. - x, beta);
		return norm * scale;
	}
}

}

Real TestSfSetParams::Params::operator()(Real x, Real z, Real Q_sq, Real ph_t_sq) const {
	Real Q = std::sqrt(Q_sq);
	Real ph_t = std::sqrt(ph_t_sq);
	Real Q_fac = std::pow(Q, -q);
	Real x_fac = L(x, alpha, beta) + f * L(x, gamma, delta);
	Real z_fac = L(z, eta, zeta);
	Real ph_t_fac = std::pow(ph_t, rho) * std::exp(-ph_t_sq / (2. * sq(sigma)));
	return N * Q_fac * x_fac * z_fac * ph_t_fac;
}

