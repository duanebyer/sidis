#ifndef SIDIS_INTERPOLATE_HPP
#define SIDIS_INTERPOLATE_HPP

#include <array>
#include <cstddef>
#include <limits>
#include <stdexcept>
#include <vector>

namespace sidis {
namespace interp {

template<typename T>
inline T linear(T f1, T f2, T x) {
	return (1. - x) * f1 + x * f2;
}

template<typename T>
inline T cubic(T f1, T f2, T f3, T f4, T x) {
	return
		2. * (1. - x) * (1. - x) * (0.5 + x) * f2
		+ 0.5 * (1. - x) * (1. - x) * x * (f3 - f1)
		+ 2. * x * x * (1.5 - x) * f3
		- 0.5 * x * x * (1. - x) * (f4 - f2);
}

template<typename T, std::size_t N>
class Grid;

/**
 * A view into a `Grid` (or another dataset).
 */
template<typename T, std::size_t N>
class GridView {
public:
	using CellIndex = std::array<std::size_t, N>;
	using Point = std::array<T, N>;

private:
	T const* _data;
	CellIndex _count;
	std::size_t _count_total;
	Point _lower;
	Point _upper;

public:
	GridView(T const* data, CellIndex count, Point lower, Point upper);

	std::size_t count(std::size_t idx) const {
		return _count[idx];
	}
	CellIndex count() const {
		return _count;
	}
	std::size_t count_total() const {
		return _count_total;
	}
	T lower(std::size_t idx) const {
		return _lower[idx];
	}
	Point lower() const {
		return _lower;
	}
	T upper(std::size_t idx) const {
		return _upper[idx];
	}
	Point upper() const {
		return _upper;
	}
	T const* data() const {
		return _data;
	}

	GridView<T, N - 1> operator[](std::size_t idx) const;
	T const& operator[](CellIndex cell_idx) const;
};

template<typename T>
class GridView<T, 0> {
public:
	using CellIndex = std::array<std::size_t, 0>;
	using Point = std::array<T, 0>;

private:
	T const* _data;

public:
	explicit GridView(T const* data) : _data(data) { }
	GridView(T const* data, CellIndex count, Point lower, Point upper) :
			_data(data) {
		static_cast<void>(count);
		static_cast<void>(lower);
		static_cast<void>(upper);
	}

	CellIndex count() const {
		return CellIndex();
	}
	std::size_t count_total() const {
		return 1;
	}
	Point lower() const {
		return Point();
	}
	Point upper() const {
		return Point();
	}
	T const* data() const {
		return _data;
	}

	T const& operator[](CellIndex cell_idx) {
		return *_data;
	}
	operator T const&() const {
		return *_data;
	}
};

/**
 * A set of N-d data on a uniform grid.
 */
template<typename T, std::size_t N>
class Grid {
public:
	using CellIndex = typename GridView<T, N>::CellIndex;
	using Point = typename GridView<T, N>::Point;

private:
	std::vector<T> _data;
	CellIndex _count;
	std::size_t _count_total;
	Point _lower;
	Point _upper;

public:
	Grid() : _data(), _count(), _count_total(0), _lower(), _upper() { }
	Grid(T const* data, CellIndex count, Point lower, Point upper);
	operator GridView<T, N>() const {
		return GridView<T, N>(_data.data(), _count, _lower, _upper);
	}

	std::size_t count(std::size_t idx) const {
		return _count[idx];
	}
	CellIndex count() const {
		return _count;
	}
	std::size_t count_total() const {
		return _count_total;
	}
	T lower(std::size_t idx) const {
		return _lower[idx];
	}
	Point lower() const {
		return _lower;
	}
	T upper(std::size_t idx) const {
		return _upper[idx];
	}
	Point upper() const {
		return _upper;
	}
	T const* data() const {
		return _data.data();
	}

	GridView<T, N - 1> operator[](std::size_t idx) const {
		return static_cast<GridView<T, N> >(*this)[idx];
	}
	T const& operator[](CellIndex cell_idx) const {
		return static_cast<GridView<T, N> >(*this)[cell_idx];
	}
};

/**
 * A function that uses linear interpolation for a regularly spaced rectangular
 * grid of N-d data.
 */
template<typename T, std::size_t N>
class LinearView {
	static_assert(
		std::is_floating_point<T>::value,
		"Cannot have a LinearView of a non-floating-point type.");

	GridView<T, N> _grid;

public:
	explicit LinearView(GridView<T, N> grid) : _grid(grid) { }
	T operator()(typename GridView<T, N>::Point x) const;
};

template<typename T>
class LinearView<T, 0> {
	static_assert(
		std::is_floating_point<T>::value,
		"Cannot have a LinearView of a non-floating-point type.");

	GridView<T, 0> _grid;

public:
	explicit LinearView(GridView<T, 0> grid) : _grid(grid) { }
	T operator()(typename GridView<T, 0>::Point x) const {
		return _grid;
	}
};

/**
 * A function that uses cubic interpolation for a regularly spaced rectangular
 * grid of N-d data.
 */
template<typename T, std::size_t N>
class CubicView {
	static_assert(
		std::is_floating_point<T>::value,
		"Cannot have a CubicView of a non-floating-point type.");

	GridView<T, N> _grid;

public:
	explicit CubicView(GridView<T, N> grid) : _grid(grid) { }
	T operator()(typename GridView<T, N>::Point x) const;
};

template<typename T>
class CubicView<T, 0> {
	static_assert(
		std::is_floating_point<T>::value,
		"Cannot have a CubicView of a non-floating-point type.");

	GridView<T, 0> _grid;

public:
	explicit CubicView(GridView<T, 0> grid) : _grid(grid) { }
	T operator()(typename GridView<T, 0>::Point x) const {
		static_cast<void>(x);
		return _grid;
	}
};

/**
 * Loads grids from an array of tuples.
 */
template<typename T, std::size_t N, std::size_t K = 1>
std::array<Grid<T, N>, K> read_grids(
	std::vector<std::array<T, N + K> > const& raw_data,
	T tolerance = 1.e2 * std::numeric_limits<T>::epsilon());

struct NotEnoughPointsError : public std::runtime_error {
	std::size_t points;
	std::size_t expected_points;
	NotEnoughPointsError(std::size_t points, std::size_t expected_points);
};

struct SingularDimensionError : public std::runtime_error {
	std::size_t dim;
	SingularDimensionError(std::size_t dim);
};

struct InvalidBoundsError : public std::runtime_error {
	InvalidBoundsError();
};

struct InvalidSpacingError : public std::runtime_error {
	InvalidSpacingError();
};

struct UnexpectedGridPointError : public std::runtime_error {
	std::size_t line_number;
	UnexpectedGridPointError(std::size_t line_number);
};

}
}

#include "sidis/extra/interpolate.ipp"

#endif

