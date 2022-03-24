#include <catch2/catch.hpp>

#include <cmath>
#include <fstream>
#include <istream>
#include <limits>
#include <sstream>
#include <utility>

#include <sidis/extra/math.hpp>
#include <sidis/extra/map.hpp>

#include "abs_matcher.hpp"
#include "rel_matcher.hpp"
#include "stream_generator.hpp"

using namespace sidis;

namespace {

template<typename T>
struct TestPair {
	T input;
	T output;
};

template<typename T>
std::istream& operator>>(std::istream& in, TestPair<T>& pair) {
	in >> pair.input;
	in >> pair.output;
	return in;
}

}

TEMPLATE_TEST_CASE(
		"√(1+x)-1 test values", "[math]",
		float, double, long double) {
	TestPair<TestType> test_pair = GENERATE(
		from_stream<TestPair<TestType>>(
			std::move(std::ifstream("data/sqrt1p_1m_vals.dat")),
			true));
	TestType x = test_pair.input;
	TestType y = test_pair.output;
	TestType y_test = math::sqrt1p_1m(x);
	std::stringstream ss;
	ss << "x = " << x;
	INFO(ss.str());
	CHECK_THAT(
		y_test,
		RelMatcher<TestType>(y, 10. * std::numeric_limits<TestType>::epsilon()));
}

TEMPLATE_TEST_CASE(
		"Dilog test values", "[math]",
		float, double, long double) {
	TestPair<TestType> test_pair = GENERATE(
		from_stream<TestPair<TestType>>(
			std::move(std::ifstream("data/dilog_vals.dat")),
			true));
	TestType x = test_pair.input;
	TestType y = test_pair.output;
	TestType y_test = math::dilog(x);
	std::stringstream ss;
	ss << "x = " << x;
	INFO(ss.str());
	CHECK_THAT(
		y_test,
		RelMatcher<TestType>(y, 10. * std::numeric_limits<TestType>::epsilon()));
}

TEST_CASE(
		"Map correctness",
		"[math]") {
	math::map::Linear map_linear;
	math::map::Inverse map_inverse(1.42);
	math::map::Log map_log(1.42);
	math::map::Decay map_decay(0.63);
	math::map::Sigmoid map_sigmoid(0.89, 0.43);
	math::map::Sigmoid2 map_sigmoid2(0.42, 0.32, 1.8, 0.45);
	Real prec = 10. * std::numeric_limits<Real>::epsilon();
	Real h = std::cbrt(std::numeric_limits<Real>::epsilon());
	Real prec_h = 1e4 * h * h;
	Real x = GENERATE(-0.64, -0.63, -0.21, 0., 0.13, 0.56, 0.93, 1.38, 1.39);
	std::stringstream ss;
	ss << "x = " << x;
	INFO(ss.str());

	// Maps should all be invertible on the interval, and calculated Jacobians
	// should be correct.
	Real x_inv, xp, xp_est, u;
	u = map_linear.u(x);
	xp_est = (map_linear.x(u + h, &xp) - map_linear.x(u - h, &xp)) / (2. * h);
	x_inv = map_linear.x(u, &xp);
	CHECK_THAT(x_inv, RelMatcher<Real>(x, prec));
	CHECK_THAT(xp, AbsMatcher<Real>(xp_est, prec_h));

	u = map_inverse.u(x);
	xp_est = (map_inverse.x(u + h, &xp) - map_inverse.x(u - h, &xp)) / (2. * h);
	x_inv = map_inverse.x(u, &xp);
	CHECK_THAT(x_inv, RelMatcher<Real>(x, prec));
	CHECK_THAT(xp, AbsMatcher<Real>(xp_est, prec_h));

	u = map_log.u(x);
	xp_est = (map_log.x(u + h, &xp) - map_log.x(u - h, &xp)) / (2. * h);
	x_inv = map_log.x(u, &xp);
	CHECK_THAT(x_inv, RelMatcher<Real>(x, prec));
	CHECK_THAT(xp, AbsMatcher<Real>(xp_est, prec_h));

	u = map_decay.u(x);
	xp_est = (map_decay.x(u + h, &xp) - map_decay.x(u - h, &xp)) / (2. * h);
	x_inv = map_decay.x(u, &xp);
	CHECK_THAT(x_inv, RelMatcher<Real>(x, prec));
	CHECK_THAT(xp, AbsMatcher<Real>(xp_est, prec_h));

	u = map_sigmoid.u(x);
	xp_est = (map_sigmoid.x(u + h, &xp) - map_sigmoid.x(u - h, &xp)) / (2. * h);
	x_inv = map_sigmoid.x(u, &xp);
	CHECK_THAT(x_inv, RelMatcher<Real>(x, prec));
	CHECK_THAT(xp, AbsMatcher<Real>(xp_est, prec_h));

	u = map_sigmoid2.u(x);
	xp_est = (map_sigmoid2.x(u + h, &xp) - map_sigmoid2.x(u - h, &xp)) / (2. * h);
	x_inv = map_sigmoid2.x(u, &xp);
	CHECK_THAT(x_inv, RelMatcher<Real>(x, prec));
	CHECK_THAT(xp, AbsMatcher<Real>(xp_est, prec_h));
}

