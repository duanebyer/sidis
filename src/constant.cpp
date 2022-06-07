#include "sidis/constant.hpp"

#include <limits>

#include "sidis/extra/interpolate.hpp"

using namespace sidis;

Real const sidis::INF = std::numeric_limits<Real>::infinity();

namespace {

Real const Q_SQ_MAX = 100.;
Real const NF = 1.;
Real const BETA_0 = 4. / 3. * NF;

Real beta(Real alpha) {
	Real alpha_r = alpha / (4. * PI);
	return 4. * PI * BETA_0 * (1. + 3. * alpha_r) * alpha_r * alpha_r;
}

interp::Grid<Real, 1> create_alpha_grid(Real Q_sq_max) {
	std::vector<Real> ln_Q_sq_vals { 2. * std::log(MASS_E) };
	std::vector<Real> alpha_vals { ALPHA_0 };
	// Generate points until a couple of points after max Q^2 is reached, to
	// ensure that cubic interpolation will work well until the end.
	std::size_t extra = 8;
	while (ln_Q_sq_vals.size() < extra
			|| ln_Q_sq_vals[ln_Q_sq_vals.size() - extra] <= std::log(Q_sq_max)) {
		// RK4 evolution.
		Real step = 0.01;
		Real ln_Q_sq_0 = ln_Q_sq_vals.back();
		Real ln_Q_sq_1 = ln_Q_sq_0 + step;
		Real alpha_0 = alpha_vals.back();
		Real alpha_p1 = beta(alpha_0);
		Real alpha_p2 = beta(alpha_0 + 0.5 * step * alpha_p1);
		Real alpha_p3 = beta(alpha_0 + 0.5 * step * alpha_p2);
		Real alpha_p4 = beta(alpha_0 + step * alpha_p3);
		Real alpha_1 = alpha_0
			+ step * (alpha_p1 + 2. * alpha_p2 + 2. * alpha_p3 + alpha_p4) / 6.;
		alpha_vals.push_back(alpha_1);
		ln_Q_sq_vals.push_back(ln_Q_sq_1);
	}
	return interp::Grid<Real, 1>(
		alpha_vals.data(),
		{ ln_Q_sq_vals.size() }, { ln_Q_sq_vals.front() }, { ln_Q_sq_vals.back() });
}

interp::Grid<Real, 1> const alpha_grid = create_alpha_grid(Q_SQ_MAX);

}

Real sidis::alpha(Real Q_sq) {
	auto cubic_view = interp::CubicView<Real, 1>(alpha_grid);
	if (Q_sq < Q_SQ_MAX) {
		// Note that if Q^2 is smaller than the electron mass, this will NaN.
		// Since this is well out of the DIS regime, this isn't a problem.
		return cubic_view({ std::log(Q_sq) });
	} else {
		// Use analytic expression to continue past the end of the grid. This is
		// the one-loop result.
		Real alpha_r_max = cubic_view({ std::log(Q_SQ_MAX) }) / (4. * PI);
		return 4. * PI * alpha_r_max
			/ (1. - alpha_r_max * BETA_0 * std::log(Q_sq / Q_SQ_MAX));
	}
}

