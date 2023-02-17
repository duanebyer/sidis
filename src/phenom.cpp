#include "sidis/phenom.hpp"

#include <cmath>

#include "sidis/constant.hpp"
#include "sidis/kinematics.hpp"
#include "sidis/extra/interpolate.hpp"

using namespace sidis;
using namespace sidis::kin;
using namespace sidis::ph;

namespace {

Real const Q_SQ_MAX = 100.;
Real const NF = 1.;
Real const BETA_0 = 4. / 3. * NF;

Real beta_qed(Real alpha) {
	Real alpha_r = alpha / (4. * PI);
	return 4. * PI * BETA_0 * (1. + 3. * alpha_r) * alpha_r * alpha_r;
}

interp::Grid<Real, 1> create_alpha_qed_grid(Real Q_sq_max) {
	std::vector<Real> ln_Q_sq_vals { 2. * std::log(MASS_E) };
	std::vector<Real> alpha_vals { ALPHA_QED_0 };
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
		Real alpha_p1 = beta_qed(alpha_0);
		Real alpha_p2 = beta_qed(alpha_0 + 0.5 * step * alpha_p1);
		Real alpha_p3 = beta_qed(alpha_0 + 0.5 * step * alpha_p2);
		Real alpha_p4 = beta_qed(alpha_0 + step * alpha_p3);
		Real alpha_1 = alpha_0
			+ step * (alpha_p1 + 2. * alpha_p2 + 2. * alpha_p3 + alpha_p4) / 6.;
		alpha_vals.push_back(alpha_1);
		ln_Q_sq_vals.push_back(ln_Q_sq_1);
	}
	return interp::Grid<Real, 1>(
		alpha_vals.data(),
		{ ln_Q_sq_vals.size() }, { ln_Q_sq_vals.front() }, { ln_Q_sq_vals.back() });
}

interp::Grid<Real, 1> const alpha_qed_grid = create_alpha_qed_grid(Q_SQ_MAX);

}

Real ph::alpha_qed(Real Q_sq) {
	auto cubic_view = interp::CubicView<Real, 1>(alpha_qed_grid);
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

Real ph::delta_vac_had(Kinematics const& kin) {
	// TODO: This vacuum hadron polarization calculation needs to be replaced
	// with a better one.
	Real alpha = 7.2973525664e-3L;
	if (kin.Q_sq < 1.) {
		return -(2.*PI)/alpha*(-1.345e-9L - 2.302e-3L*std::log(1. + 4.091L*kin.Q_sq));
	} else if (kin.Q_sq < 64.) {
		return -(2.*PI)/alpha*(-1.512e-3L - 2.822e-3L*std::log(1. + 1.218L*kin.Q_sq));
	} else {
		return -(2.*PI)/alpha*(-1.1344e-3L - 3.0680e-3L*std::log(1. + 0.99992L*kin.Q_sq));
	}
}

Phenom::Phenom(Kinematics const& kin) :
	alpha_qed(ph::alpha_qed(kin.Q_sq)),
	delta_vac_had(ph::delta_vac_had(kin)) { }
Phenom::Phenom(Real alpha_qed, Kinematics const& kin) :
	alpha_qed(alpha_qed),
	delta_vac_had(ph::delta_vac_had(kin)) { }
Phenom::Phenom(Real alpha_qed, Real delta_vac_had) :
	alpha_qed(alpha_qed),
	delta_vac_had(delta_vac_had) { }

