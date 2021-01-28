#ifndef SIDIS_INTERPOLATE_IPP
#define SIDIS_INTERPOLATE_IPP

#include "sidis/extra/interpolate.hpp"

#include <cmath>
#include <cstddef>
#include <limits>
#include <stdexcept>
#include <string>
#include <type_traits>

namespace sidis {
namespace interp {

template<typename T, std::size_t N>
GridView<T, N>::GridView(
		T const* data,
		CellIndex count,
		Point lower,
		Point upper) :
		_data(data),
		_count(count),
		_lower(lower),
		_upper(upper) {
	_count_total = 1;
	for (std::size_t idx = 0; idx < N; ++idx) {
		_count_total *= _count[idx];
		if (!std::isfinite(lower[idx]) || !std::isfinite(upper[idx])) {
			throw InvalidBoundsError();
		}
		if (lower[idx] > upper[idx]) {
			throw InvalidBoundsError();
		}
		if (_count[idx] <= 1) {
			throw SingularDimensionError(idx);
		}
	}
}

template<typename T, std::size_t N>
GridView<T, N - 1> GridView<T, N>::operator[](std::size_t idx) const {
	std::size_t width = _count_total / _count[0];
	typename GridView<T, N - 1>::CellIndex next_count;
	typename GridView<T, N - 1>::Point next_lower;
	typename GridView<T, N - 1>::Point next_upper;
	for (std::size_t idx = 1; idx < N; ++idx) {
		next_count[idx - 1] = _count[idx];
		next_lower[idx - 1] = _lower[idx];
		next_upper[idx - 1] = _upper[idx];
	}
	return GridView<T, N - 1>(
		_data + idx * width,
		next_count,
		next_lower,
		next_upper);
}

template<typename T, std::size_t N>
T const& GridView<T, N>::operator[](CellIndex cell_idx) const {
	typename GridView<T, N - 1>::CellIndex next_cell_idx;
	for (std::size_t idx = 1; idx < N; ++idx) {
		next_cell_idx[idx - 1] = cell_idx[idx];
	}
	return (*this)[cell_idx[0]][next_cell_idx];
}

template<typename T, std::size_t N>
Grid<T, N>::Grid(T const* data, CellIndex count, Point lower, Point upper) :
		_count(count),
		_lower(lower),
		_upper(upper) {
	_count_total = 1;
	for (std::size_t idx = 0; idx < N; ++idx) {
		_count_total *= _count[idx];
		if (!std::isfinite(lower[idx]) || !std::isfinite(upper[idx])) {
			throw InvalidBoundsError();
		}
		if (lower[idx] > upper[idx]) {
			throw InvalidBoundsError();
		}
		if (_count[idx] <= 1) {
			throw SingularDimensionError(idx);
		}
	}
	_data = std::vector<T>(data, data + _count_total);
}

template<typename T, std::size_t N>
T LinearView<T, N>::operator()(typename GridView<T, N>::Point x) const {
	T x_diff = x[0] - _grid.lower(0);
	T width = _grid.upper(0) - _grid.lower(0);
	if (std::isnan(x_diff)) {
		return x_diff;
	} else if (x_diff < 0. || x_diff >= width) {
		return std::numeric_limits<T>::quiet_NaN();
	}
	T idx_float;
	T x_rel = std::modf((_grid.count(0) - 1) * (x_diff / width), &idx_float);
	std::size_t idx = (std::size_t) idx_float;
	typename GridView<T, N - 1>::Point x_next;
	for (std::size_t idx = 1; idx < N; ++idx) {
		x_next[idx - 1] = x[idx];
	}
	T linear_1 = LinearView<T, N - 1>(_grid[idx])(x_next);
	T linear_2 = LinearView<T, N - 1>(_grid[idx + 1])(x_next);
	return linear(linear_1, linear_2, x_rel);
}

template<typename T, std::size_t N>
T CubicView<T, N>::operator()(typename GridView<T, N>::Point x) const {
	T x_diff = x[0] - _grid.lower(0);
	T width = _grid.upper(0) - _grid.lower(0);
	if (std::isnan(x_diff)) {
		return x_diff;
	} else if (x_diff < 0. || x_diff >= width) {
		return std::numeric_limits<T>::quiet_NaN();
	}
	T idx_float;
	T x_rel = std::modf((_grid.count(0) - 1) * (x_diff / width), &idx_float);
	std::size_t idx = (std::size_t) idx_float;
	typename GridView<T, N - 1>::Point x_next;
	for (std::size_t idx = 1; idx < N; ++idx) {
		x_next[idx - 1] = x[idx];
	}
	T cubic_0 = idx == 0 ? 0. :
		CubicView<T, N - 1>(_grid[idx - 1])(x_next);
	T cubic_1 = CubicView<T, N - 1>(_grid[idx])(x_next);
	T cubic_2 = CubicView<T, N - 1>(_grid[idx + 1])(x_next);
	T cubic_3 = idx >= _grid.count(0) - 2 ? 0. :
		CubicView<T, N - 1>(_grid[idx + 2])(x_next);
	return cubic(cubic_0, cubic_1, cubic_2, cubic_3, x_rel);
}


template<typename T, std::size_t N, std::size_t K>
inline std::array<Grid<T, N>, K> read_grids(
		std::vector<std::array<T, N + K> > const& raw_data,
		T tolerance) {
	std::vector<typename GridView<T, N>::Point> grid_points(raw_data.size());
	std::vector<std::array<T, K> > data(raw_data.size());
	for (std::size_t row_idx = 0; row_idx < raw_data.size(); ++row_idx) {
		for (std::size_t dim_idx = 0; dim_idx < N; ++dim_idx) {
			grid_points[row_idx][dim_idx] = raw_data[row_idx][dim_idx];
		}
		for (std::size_t col_idx = 0; col_idx < K; ++col_idx) {
			data[row_idx][col_idx] = raw_data[row_idx][N + col_idx];
		}
	}

	// If there are not enough grid points to form one unit cell, then give up.
	if (grid_points.size() < (1 << N)) {
		throw NotEnoughPointsError(grid_points.size(), 1 << N);
	}

	std::size_t count_total = 1;
	std::array<bool, N> filled_cols = { false };
	std::array<std::size_t, N> dim_permute_map;
	typename GridView<T, N>::CellIndex counts;
	typename GridView<T, N>::Point lower;
	typename GridView<T, N>::Point upper;
	for (std::size_t dim_idx = 0; dim_idx < N; ++dim_idx) {
		std::size_t sub_count = grid_points.size() / count_total;
		if (sub_count * count_total != grid_points.size()) {
			throw NotEnoughPointsError(
				grid_points.size(),
				(sub_count + 1) * count_total);
		}
		// Find the column that steps by `count_total`.
		std::size_t step_dim_idx = 0;
		while (true) {
			if (step_dim_idx == N) {
				throw NotEnoughPointsError(
					grid_points.size(),
					2 * grid_points.size());
			}
			if (filled_cols[step_dim_idx]) {
				step_dim_idx += 1;
				continue;
			}
			T initial_point = grid_points[0][step_dim_idx];
			T final_point = grid_points[count_total][step_dim_idx];
			if (initial_point == final_point) {
				step_dim_idx += 1;
				continue;
			}
			bool same = true;
			for (std::size_t row_idx = 0; row_idx < count_total; ++row_idx) {
				if (grid_points[row_idx][step_dim_idx] != initial_point) {
					same = false;
					break;
				}
			}
			if (!same) {
				step_dim_idx += 1;
				continue;
			}
			filled_cols[step_dim_idx] = true;
			dim_permute_map[dim_idx] = step_dim_idx;
			break;
		}

		// Use the column to fill a vector of the grid spacings.
		std::vector<T> grid_planes;
		for (std::size_t group_idx = 0; group_idx < sub_count; ++group_idx) {
			T next = grid_points[group_idx * count_total][step_dim_idx];
			if (grid_planes.empty() || next != grid_planes[0]) {
				grid_planes.push_back(next);
			} else {
				break;
			}
		}
		std::size_t next_count = grid_planes.size();
		counts[step_dim_idx] = next_count;
		if (next_count <= 1) {
			throw SingularDimensionError(step_dim_idx);
		}

		// Check that the data in the column is valid.
		for (std::size_t row_idx = 0; row_idx < grid_points.size(); ++row_idx) {
			T next = grid_points[row_idx][step_dim_idx];
			std::size_t group_idx = row_idx / count_total;
			if (next != grid_planes[group_idx % next_count]) {
				throw UnexpectedGridPointError(row_idx);
			}
		}

		// Check that the lower bound is less than the upper bound.
		T next_lower = grid_planes.front();
		T next_upper = grid_planes.back();
		lower[step_dim_idx] = next_lower;
		upper[step_dim_idx] = next_upper;
		if (!std::isfinite(next_lower) || !std::isfinite(next_upper)) {
			throw InvalidBoundsError();
		}
		if (next_lower >= next_upper) {
			throw InvalidBoundsError();
		}

		// Check that the grids are spaced approximately evenly.
		T spacing = (next_upper - next_lower) / (next_count - 1);
		for (std::size_t pl_idx = 0; pl_idx < next_count; ++pl_idx) {
			T plane = grid_planes[pl_idx];
			T uniform_plane = next_lower + pl_idx * spacing;
			T rel_tol = tolerance;
			if (std::abs(plane - uniform_plane) > std::abs(rel_tol * spacing)) {
				throw InvalidSpacingError();
			}
		}

		// Increase count for next round.
		count_total *= next_count;
	}

	// Use the grid information to reorder the data.
	std::array<std::vector<T>, K> data_transposed;
	std::array<std::size_t, N> count_totals = { count_total / counts[0] };
	for (std::size_t dim_idx = 1; dim_idx < N; ++dim_idx) {
		count_totals[dim_idx] = count_totals[dim_idx - 1] / counts[dim_idx];
	}
	for (std::size_t col_idx = 0; col_idx < K; ++col_idx) {
		data_transposed[col_idx].resize(data.size());
	}
	for (std::size_t row_idx = 0; row_idx < data.size(); ++row_idx) {
		std::size_t new_row_idx = 0;
		std::size_t new_count = 1;
		for (std::size_t dim_idx = 0; dim_idx < N; ++dim_idx) {
			std::size_t new_dim_idx = dim_permute_map[dim_idx];
			std::size_t rel_idx = (row_idx / new_count) % counts[new_dim_idx];
			std::size_t old_count = count_totals[new_dim_idx];
			new_row_idx += rel_idx * old_count;
			new_count *= counts[new_dim_idx];
		}
		for (std::size_t col_idx = 0; col_idx < K; ++col_idx) {
			data_transposed[col_idx][new_row_idx] = data[row_idx][col_idx];
		}
	}

	std::array<Grid<T, N>, K> grids;
	for (std::size_t col_idx = 0; col_idx < K; ++col_idx) {
		grids[col_idx] = Grid<T, N>(
			data_transposed[col_idx].data(),
			counts, lower, upper);
	}

	return grids;
}

inline NotEnoughPointsError::NotEnoughPointsError(
	std::size_t points,
	std::size_t expected_points) :
	std::runtime_error(
		"Not enough points to construct the grid (found "
		+ std::to_string(points) + ", expected "
		+ std::to_string(expected_points) + ")"),
	points(points),
	expected_points(expected_points) { }

inline SingularDimensionError::SingularDimensionError(std::size_t dim) :
	std::runtime_error("Grid is singular in dimension " + std::to_string(dim)),
	dim(dim) { }

inline InvalidBoundsError::InvalidBoundsError() :
	std::runtime_error("Lower grid bound must be smaller than upper bound") { }

inline InvalidSpacingError::InvalidSpacingError() :
	std::runtime_error("Grid must be spaced uniformly") { }

inline UnexpectedGridPointError::UnexpectedGridPointError(
	std::size_t line_number) :
	std::runtime_error(
		"Unexpected grid point at line " + std::to_string(line_number)),
	line_number(line_number) { }

}
}

#endif

