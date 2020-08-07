#include <catch2/catch.hpp>

#include <cmath>
#include <tuple>

#include <sidis/math.hpp>

using namespace sidis;

TEMPLATE_TEST_CASE(
		"Dilog test values", "[math]",
		float, double, long double) {
	using pair = std::tuple<TestType, TestType>;
	auto target = GENERATE(
		pair {   0.00,  0                     },
		pair {  -0.20, -0.1908                },
		pair {   0.25,  0.267653              },
		pair {  -0.35, -0.32337               },
		pair {   0.40,  0.449283              },
		pair {  -0.50, -0.448414              },
		pair {   0.80,  1.07479               },
		pair {  -0.90, -0.752163              },
		pair {   0.99,  1.58863               },
		pair {   1.00,  1.6449340668482264365 },
		pair {   1.01,  1.70073               },
		pair {  -1.50, -1.14738               },
		pair {   2.00,  2.4674011002723396547 },
		pair {  -4.00, -2.3699397969983658320 },
		pair {   5.00,  1.7837191612666306277 },
		pair { -16.00, -5.4270085310681956913 },
		pair { 100.00, -7.3239531990004822427 });

	TestType x = std::get<0>(target);
	TestType y = std::get<1>(target);
	TestType y_test = math::dilog(x);
	INFO("x = " + std::to_string(x));
	CHECK(y_test == Approx(y));
}

