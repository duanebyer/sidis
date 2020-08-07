#include <catch2/catch.hpp>

#include <cmath>
#include <tuple>

#include <sidis/math.hpp>

using namespace sidis;

TEMPLATE_TEST_CASE(
		"Dilog test values", "[math]",
		float, double, long double) {
	using Pair = std::tuple<TestType, TestType>;
	auto target = GENERATE(
		Pair {   0.00,  0                     },
		Pair {  -0.20, -0.1908                },
		Pair {   0.25,  0.267653              },
		Pair {  -0.35, -0.32337               },
		Pair {   0.40,  0.449283              },
		Pair {  -0.50, -0.448414              },
		Pair {   0.80,  1.07479               },
		Pair {  -0.90, -0.752163              },
		Pair {   0.99,  1.58863               },
		Pair {   1.00,  1.6449340668482264365 },
		Pair {   1.01,  1.70073               },
		Pair {  -1.50, -1.14738               },
		Pair {   2.00,  2.4674011002723396547 },
		Pair {  -4.00, -2.3699397969983658320 },
		Pair {   5.00,  1.7837191612666306277 },
		Pair { -16.00, -5.4270085310681956913 },
		Pair { 100.00, -7.3239531990004822427 });

	TestType x = std::get<0>(target);
	TestType y = std::get<1>(target);
	TestType y_test = math::dilog(x);
	INFO("x = " + std::to_string(x));
	CHECK(y_test == Approx(y));
}

