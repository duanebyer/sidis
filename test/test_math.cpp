#include <catch2/catch.hpp>

#include <cmath>
#include <fstream>
#include <istream>
#include <limits>
#include <sstream>
#include <utility>

#include <sidis/math.hpp>

#include "file_generator.hpp"
#include "rel_matcher.hpp"

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

