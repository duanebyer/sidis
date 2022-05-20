#ifndef SIDIS_INTEGRATE_HPP
#define SIDIS_INTEGRATE_HPP

#include <algorithm>
#include <array>
#include <memory>
#include <new>

#include <cubature.hpp>

#include <gsl/gsl_monte.h>
#include <gsl/gsl_rng.h>
#include <gsl/gsl_monte_miser.h>
#include <gsl/gsl_monte_plain.h>
#include <gsl/gsl_monte_vegas.h>

#include "sidis/integ_params.hpp"
#include "sidis/numeric.hpp"

namespace sidis {
namespace math {

/// Integrates over a hyper-cube.
template<std::size_t D, typename F>
EstErr integrate(F func,
		std::array<Real, D> lower,
		std::array<Real, D> upper,
		IntegParams params) {
	if (params.method == IntegMethod::CUBATURE) {
		cubature::EstErr<Real> result = cubature::cubature<D>(
			func,
			lower, upper,
			params.num_evals,
			params.err_abs,
			params.err_rel);
		return { result.val, result.err };
	} else if (
			params.method == IntegMethod::MC_PLAIN
			|| params.method == IntegMethod::MISER
			|| params.method == IntegMethod::VEGAS) {
		EstErr result;
		// TODO: How to deal with the case where `Real` is not `double`? This
		// code will fail to compile as written.
		auto func_ptr = [](Real* xs, std::size_t, void* params) {
			F* func = static_cast<F*>(params);
			std::array<Real, D> point;
			std::move(xs, xs + D, point.begin());
			return (*func)(point);
		};
		gsl_monte_function func_gsl { func_ptr, D, static_cast<void*>(&func) };
		// TODO: Which random number generator should be used by default?
		std::unique_ptr<gsl_rng, void(*)(gsl_rng*)> rnd(
			gsl_rng_alloc(gsl_rng_taus2), gsl_rng_free);
		if (rnd == nullptr) {
			throw std::bad_alloc();
		}
		if (params.method == IntegMethod::MC_PLAIN) {
			std::unique_ptr<
					gsl_monte_plain_state,
					void(*)(gsl_monte_plain_state*)> state(
				gsl_monte_plain_alloc(D), &gsl_monte_plain_free);
			if (state == nullptr) {
				throw std::bad_alloc();
			}
			gsl_monte_plain_integrate(
				&func_gsl,
				lower.data(),
				upper.data(),
				D,
				params.num_evals,
				rnd.get(),
				state.get(),
				&result.val,
				&result.err);
		} else if (params.method == IntegMethod::MISER) {
			std::unique_ptr<
					gsl_monte_miser_state,
					void(*)(gsl_monte_miser_state*)> state(
				gsl_monte_miser_alloc(D), &gsl_monte_miser_free);
			if (state == nullptr) {
				throw std::bad_alloc();
			}
			gsl_monte_miser_integrate(
				&func_gsl,
				lower.data(),
				upper.data(),
				D,
				params.num_evals,
				rnd.get(),
				state.get(),
				&result.val,
				&result.err);
		} else if (params.method == IntegMethod::VEGAS) {
			std::unique_ptr<
					gsl_monte_vegas_state,
					void(*)(gsl_monte_vegas_state*)> state(
				gsl_monte_vegas_alloc(D), &gsl_monte_vegas_free);
			if (state == nullptr) {
				throw std::bad_alloc();
			}
			// TODO: Divide into warm-up and actual runs. Might be worth
			// gradually increasing the number of runs until the error is below
			// the threshold.
			gsl_monte_vegas_integrate(
				&func_gsl,
				lower.data(),
				upper.data(),
				D,
				params.num_evals,
				rnd.get(),
				state.get(),
				&result.val,
				&result.err);
		}
		return result;
	} else {
		return { 0., 0. };
	}
}

/// Integrates over an interval.
template<typename F>
EstErr integrate(F func, Real lower, Real upper, IntegParams params) {
	auto func_arr = [&](std::array<Real, 1> x) { return func(x[0]); };
	return integrate<1>(func_arr, { lower }, { upper }, params);
}

template<typename F>
Real riemann(F f, Real a, Real b, unsigned n) {
	Real delta = (b - a) / n;
	Real result = 0.;
	for (unsigned i = 0; i < n; ++i) {
		Real x = (i * b + (n - i) * a + 0.5 * (b - a)) / n;
		result += f(x);
	}
	return result * delta;
}

template<typename F>
Real trapezoid(F f, Real a, Real b, unsigned n) {
	Real delta = (b - a) / n;
	Real result = 0.5 * (f(a) + f(b));
	for (unsigned i = 1; i <= n - 1; ++i) {
		Real x = (i * b + (n - i) * a) / n;
		result += f(x);
	}
	return result * delta;
}

}
}

#endif

