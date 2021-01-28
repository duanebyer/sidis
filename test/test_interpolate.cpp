#include <catch2/catch.hpp>

#include <array>
#include <cmath>
#include <fstream>
#include <limits>
#include <sstream>
#include <vector>

#include <sidis/numeric.hpp>
#include <sidis/extra/interpolate.hpp>

#include "rel_matcher.hpp"

using namespace sidis;

namespace {

template<typename T, std::size_t N>
std::istream& operator>>(
		std::istream& in,
		std::vector<std::array<T, N> >& data) {
	data.clear();
	while (in) {
		std::array<T, N> next;
		for (std::size_t idx = 0; idx < N; ++idx) {
			in >> next[idx];
		}
		if (in) {
			data.push_back(next);
		}
	}
	return in;
}

}

TEST_CASE(
		"Grid reading from array test",
		"[interp]") {
	std::ifstream data_file("data/grid_1_vals.dat");
	std::vector<std::array<double, 3 + 2> > data;
	data_file >> data;
	std::array<interp::Grid<double, 3>, 2> grids
		= interp::read_grids<double, 3, 2>(data);
	interp::Grid<double, 3> grid_1 = grids[0];
	interp::Grid<double, 3> grid_2 = grids[1];

	// Check sizes of the grids.
	CHECK(grid_1.count() == std::array<std::size_t, 3>{ 2, 3, 4 });
	CHECK(grid_2.count() == std::array<std::size_t, 3>{ 2, 3, 4 });
	CHECK(grid_1.count_total() == 24);
	CHECK(grid_2.count_total() == 24);

	// Check bounds of the grids.
	CHECK(grid_1.lower() == std::array<double, 3>{ 0., 0., 0. });
	CHECK(grid_1.upper() == std::array<double, 3>{ 1., 2., 3. });
	CHECK(grid_2.lower() == std::array<double, 3>{ 0., 0., 0. });
	CHECK(grid_2.upper() == std::array<double, 3>{ 1., 2., 3. });

	// Check all eight corners of the grids.
	CHECK(grid_1[{ 0, 0, 0 }] == -0.8);
	CHECK(grid_1[{ 1, 0, 0 }] == -2.3);
	CHECK(grid_1[{ 0, 2, 0 }] == -1.2);
	CHECK(grid_1[{ 0, 0, 3 }] == -2.0);
	CHECK(grid_1[{ 0, 2, 3 }] == 0.6);
	CHECK(grid_1[{ 1, 0, 3 }] == 0.3);
	CHECK(grid_1[{ 1, 2, 0 }] == 0.7);
	CHECK(grid_1[{ 1, 2, 3 }] == -0.1);

	CHECK(grid_2[{ 0, 0, 0 }] == -0.3);
	CHECK(grid_2[{ 1, 0, 0 }] == 0.3);
	CHECK(grid_2[{ 0, 2, 0 }] == 1.3);
	CHECK(grid_2[{ 0, 0, 3 }] == -0.6);
	CHECK(grid_2[{ 0, 2, 3 }] == -1.6);
	CHECK(grid_2[{ 1, 0, 3 }] == -1.4);
	CHECK(grid_2[{ 1, 2, 0 }] == 0.4);
	CHECK(grid_2[{ 1, 2, 3 }] == -0.2);

	// Check some random points inside the grids.
	CHECK(grid_1[{ 1, 1, 2 }] == -0.5);
	CHECK(grid_1[{ 0, 2, 1 }] == 0.3);
	CHECK(grid_1[{ 0, 1, 3 }] == -0.8);

	CHECK(grid_2[{ 1, 1, 2 }] == 0.5);
	CHECK(grid_2[{ 0, 2, 1 }] == -2.1);
	CHECK(grid_2[{ 0, 1, 3 }] == 2.4);

	// Check reading using subgrids.
	CHECK(grid_1[1][1][2] == -0.5);
	CHECK(grid_1[0][2][1] == 0.3);
	CHECK(grid_1[0][1][3] == -0.8);

	CHECK(grid_2[1][1][2] == 0.5);
	CHECK(grid_2[0][2][1] == -2.1);
	CHECK(grid_2[0][1][3] == 2.4);
}

TEST_CASE(
		"Grid reading from single cell test",
		"[interp]") {
	std::ifstream data_file("data/grid_2_vals.dat");
	std::vector<std::array<double, 3 + 1> > data;
	data_file >> data;
	interp::Grid<double, 3> grid = interp::read_grids<double, 3, 1>(data)[0];

	// Check sizes.
	CHECK(grid.count() == std::array<std::size_t, 3>{ 2, 2, 2 });
	CHECK(grid.count_total() == 8);

	// Check the corners.
	CHECK(grid[{ 0, 0, 0 }] == 1.);
	CHECK(grid[{ 0, 0, 1 }] == 2.);
	CHECK(grid[{ 0, 1, 0 }] == 3.);
	CHECK(grid[{ 0, 1, 1 }] == 4.);
	CHECK(grid[{ 1, 0, 0 }] == 5.);
	CHECK(grid[{ 1, 0, 1 }] == 6.);
	CHECK(grid[{ 1, 1, 0 }] == 7.);
	CHECK(grid[{ 1, 1, 1 }] == 8.);
}

TEST_CASE(
		"Grid reading verification tests",
		"[interp]") {
	// Try loading some grids that have issues with them, and make sure that the
	// problems are caught correctly.
	std::ifstream data_bad_1_file("data/grid_bad_1_vals.dat");
	std::ifstream data_bad_2_file("data/grid_bad_2_vals.dat");
	std::ifstream data_bad_3_file("data/grid_bad_3_vals.dat");
	std::ifstream data_bad_4_file("data/grid_bad_4_vals.dat");
	std::vector<std::array<double, 3> > data_bad_1;
	std::vector<std::array<double, 3> > data_bad_2;
	std::vector<std::array<double, 3> > data_bad_3;
	std::vector<std::array<double, 3> > data_bad_4;
	data_bad_1_file >> data_bad_1;
	data_bad_2_file >> data_bad_2;
	data_bad_3_file >> data_bad_3;
	data_bad_4_file >> data_bad_4;

	CHECK_THROWS_AS(
		(interp::read_grids<double, 3, 0>)(data_bad_1),
		interp::NotEnoughPointsError);
	CHECK_THROWS_AS(
		(interp::read_grids<double, 3, 0>)(data_bad_2),
		interp::InvalidSpacingError);
	CHECK_THROWS_AS(
		(interp::read_grids<double, 3, 0>)(data_bad_3),
		interp::UnexpectedGridPointError);
	CHECK_THROWS_AS(
		(interp::read_grids<double, 3, 0>)(data_bad_4),
		interp::UnexpectedGridPointError);
}

TEST_CASE(
		"Linear interpolation on grids tests",
		"[interp]") {
	std::ifstream data_file("data/grid_3_vals.dat");
	std::vector<std::array<double, 3 + 2> > data;
	data_file >> data;
	std::array<interp::Grid<double, 3>, 2> grids
		= interp::read_grids<double, 3, 2>(data);
	interp::LinearView<double, 3> linear(grids[0]);
	interp::LinearView<double, 3> linear_const(grids[1]);

	RelMatcher<double> matcher_const(
		2.,
		1e2 * std::numeric_limits<double>::epsilon());
	// Check in bounds.
	CHECK_THAT(linear_const({ 12., 0.8, 2.5 }), matcher_const);
	CHECK_THAT(linear_const({ 18., 1.2, 7.5 }), matcher_const);
	CHECK_THAT(linear_const({ 24., 2.0, 2.1 }), matcher_const);
	CHECK_THAT(linear_const({ 10., 0.0, 2.0 }), matcher_const);
	CHECK_THAT(linear_const({ 20., 1.0, 4.0 }), matcher_const);
	CHECK_THAT(linear_const({ 28., 1.1, 4.9 }), matcher_const);
	CHECK_THAT(linear_const({ 29., 2.9, 7.9 }), matcher_const);
	// Check out of bounds.
	CHECK(std::isnan(linear_const({ 31.0,  2.0, 4.0 })));
	CHECK(std::isnan(linear_const({  9.9,  2.0, 4.0 })));
	CHECK(std::isnan(linear_const({ 20.0,  3.1, 4.0 })));
	CHECK(std::isnan(linear_const({ 20.0, -0.1, 4.0 })));
	CHECK(std::isnan(linear_const({ 20.0,  2.0, 8.1 })));
	CHECK(std::isnan(linear_const({ 20.0,  2.0, 1.9 })));
}

TEST_CASE(
		"Cubic interpolation on grids tests",
		"[interp]") {
	std::ifstream data_file("data/grid_3_vals.dat");
	std::vector<std::array<double, 3 + 2> > data;
	data_file >> data;
	std::array<interp::Grid<double, 3>, 2> grids
		= interp::read_grids<double, 3, 2>(data);
	interp::CubicView<double, 3> cubic(grids[0]);
	interp::CubicView<double, 3> cubic_const(grids[1]);

	// The tolerance on the matcher is high because near the edges, the cubic
	// interpolator is not very constant.
	RelMatcher<double> matcher_const(
		2.,
		0.5);
	// Check in bounds.
	CHECK_THAT(cubic_const({ 12., 0.8, 2.5 }), matcher_const);
	CHECK_THAT(cubic_const({ 18., 1.2, 7.5 }), matcher_const);
	CHECK_THAT(cubic_const({ 24., 2.0, 2.1 }), matcher_const);
	CHECK_THAT(cubic_const({ 10., 0.0, 2.0 }), matcher_const);
	CHECK_THAT(cubic_const({ 20., 1.0, 4.0 }), matcher_const);
	CHECK_THAT(cubic_const({ 28., 1.1, 4.9 }), matcher_const);
	CHECK_THAT(cubic_const({ 29., 2.9, 7.9 }), matcher_const);
	// Check out of bounds.
	CHECK(std::isnan(cubic_const({ 31.0,  2.0, 4.0 })));
	CHECK(std::isnan(cubic_const({  9.9,  2.0, 4.0 })));
	CHECK(std::isnan(cubic_const({ 20.0,  3.1, 4.0 })));
	CHECK(std::isnan(cubic_const({ 20.0, -0.1, 4.0 })));
	CHECK(std::isnan(cubic_const({ 20.0,  2.0, 8.1 })));
	CHECK(std::isnan(cubic_const({ 20.0,  2.0, 1.9 })));
}

double test_function(double x, double y) {
	double pi = 3.1415926;
	return std::sin(2. * pi * x) * std::cos(2. * pi * y);
}

TEST_CASE(
		"Interpolation to approximate function tests",
		"[interp]") {
	std::vector<double> data;
	std::size_t count = 101;
	double length = 1.;
	for (std::size_t x_idx = 0; x_idx < count; ++x_idx) {
		for (std::size_t y_idx = 0; y_idx < count; ++y_idx) {
			double x = length * x_idx / (count - 1);
			double y = length * y_idx / (count - 1);
			data.push_back(test_function(x, y));
		}
	}
	interp::Grid<double, 2> grid(
		data.data(),
		{ count, count },
		{ 0., 0. },
		{ length, length });
	interp::LinearView<double, 2> linear(grid);
	interp::CubicView<double, 2> cubic(grid);

	// Compare to the test function at various points.
	using Pair = std::array<double, 2>;
	Pair pair = GENERATE(
		Pair{ 0.1283, 0.2853 },
		Pair{ 0.2944, 0.2944 },
		Pair{ 0.2043, 0.1012 },
		Pair{ 0.5000, 0.5000 },
		Pair{ 0.1103, 0.9821 },
		Pair{ 0.9384, 0.1839 });
	RelMatcher<double> matcher(
		test_function(pair[0], pair[1]),
		1e1 * length / count);
	std::stringstream ss;
	ss << "x = " << pair[0] << ", y = " << pair[1];
	INFO(ss.str());
	CHECK_THAT(linear(pair), matcher);
	CHECK_THAT(cubic(pair), matcher);
}

