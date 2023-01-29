#ifndef SIDIS_INTERPOLATE_HPP
#define SIDIS_INTERPOLATE_HPP

#include <array>
#include <cstddef>
#include <limits>
#include <stdexcept>
#include <vector>

namespace sidis {
namespace interp {

/**
 * \defgroup InterpGroup Interpolation
 * Types and methods for interpolation of gridded data.
 */
/// \{

/// Linearly interpolate between two endpoints, with \p x between 0 and 1.
template<typename T>
inline T linear(T f1, T f2, T x) {
	return (1. - x) * f1 + x * f2;
}

/// Cubic interpolation using \p f2 and \p f3 as the endpoints, and \p f1 and
/// \p f4 as the anchor points beyond the endpoints. \p x between 0 and 1.
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
 * A view into a Grid (or another dataset). The difference between a GridView
 * and a Grid is that a Grid owns its data, while a GridView provides a
 * Grid-like interface into an already existing set of data.
 */
template<typename T, std::size_t N>
class GridView final {
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
	/// Construct a GridView out of \p data. The GridView spans the hyper-cube
	/// from \p lower to \p upper. The number of data points in each dimension
	/// is determined by the \p count array.
	GridView(T const* data, CellIndex count, Point lower, Point upper);

	/// Number of data points in dimension \p idx.
	std::size_t count(std::size_t idx) const {
		return _count[idx];
	}
	/// Number of data points in each dimension.
	CellIndex count() const {
		return _count;
	}
	/// Total number of data points.
	std::size_t count_total() const {
		return _count_total;
	}
	/// Lower bound on the hyper-cube spanned by the data.
	T lower(std::size_t idx) const {
		return _lower[idx];
	}
	/// \copydoc GridView::lower()
	Point lower() const {
		return _lower;
	}
	/// Upper bound on the hyper-cube spanned by the data.
	T upper(std::size_t idx) const {
		return _upper[idx];
	}
	/// \copydoc GridView::upper()
	Point upper() const {
		return _upper;
	}
	/// Access to the contigious underlying data.
	T const* data() const {
		return _data;
	}

	/// Access a slice of the data.
	GridView<T, N - 1> operator[](std::size_t idx) const;
	/// Access an element from the grid.
	T const& operator[](CellIndex cell_idx) const;
};

/// Base specialization of GridView.
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
		static_cast<void>(cell_idx);
		return *_data;
	}
	operator T const&() const {
		return *_data;
	}
};

/**
 * A set of N-d data on a uniform grid. A Grid is an owned version of a
 * GridView.
 *
 * \sa GridView
 */
template<typename T, std::size_t N>
class Grid final {
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

	/// \copydoc GridView::count()
	std::size_t count(std::size_t idx) const {
		return _count[idx];
	}
	/// \copydoc GridView::count()
	CellIndex count() const {
		return _count;
	}
	/// \copydoc GridView::count_total()
	std::size_t count_total() const {
		return _count_total;
	}
	/// \copydoc GridView::lower()
	T lower(std::size_t idx) const {
		return _lower[idx];
	}
	/// \copydoc GridView::lower()
	Point lower() const {
		return _lower;
	}
	/// \copydoc GridView::upper()
	T upper(std::size_t idx) const {
		return _upper[idx];
	}
	/// \copydoc GridView::upper()
	Point upper() const {
		return _upper;
	}
	/// \copydoc GridView::data()
	T const* data() const {
		return _data.data();
	}

	/// \copydoc GridView::operator[]()
	GridView<T, N - 1> operator[](std::size_t idx) const {
		return static_cast<GridView<T, N> >(*this)[idx];
	}
	/// \copydoc GridView::operator[]()
	T const& operator[](CellIndex cell_idx) const {
		return static_cast<GridView<T, N> >(*this)[cell_idx];
	}
};

/**
 * A function that uses linear interpolation for a regularly spaced rectangular
 * grid of N-d data.
 */
template<typename T, std::size_t N>
class LinearView final {
	static_assert(
		std::is_floating_point<T>::value,
		"Cannot have a LinearView of a non-floating-point type.");

	GridView<T, N> _grid;

public:
	explicit LinearView(GridView<T, N> grid) : _grid(grid) { }
	/// Interpolate from the underlying GridView.
	T operator()(typename GridView<T, N>::Point x) const;
};

/// Base specialization of LinearView.
template<typename T>
class LinearView<T, 0> final {
	static_assert(
		std::is_floating_point<T>::value,
		"Cannot have a LinearView of a non-floating-point type.");

	GridView<T, 0> _grid;

public:
	explicit LinearView(GridView<T, 0> grid) : _grid(grid) { }
	T operator()(typename GridView<T, 0>::Point x) const {
		static_cast<void>(x);
		return _grid;
	}
};

/**
 * A function that uses cubic interpolation for a regularly spaced rectangular
 * grid of N-d data.
 */
template<typename T, std::size_t N>
class CubicView final {
	static_assert(
		std::is_floating_point<T>::value,
		"Cannot have a CubicView of a non-floating-point type.");

	GridView<T, N> _grid;

public:
	explicit CubicView(GridView<T, N> grid) : _grid(grid) { }
	/// Interpolate from the underlying GridView.
	T operator()(typename GridView<T, N>::Point x) const;
};

/// Base specialization of CubicView.
template<typename T>
class CubicView<T, 0> final {
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

/// Loads grids from an array of tuples. The provided data points must be
/// provided either in row-major order or column-major order. For each data
/// point, the first \p N numbers give the coordinates, and the next \p K
/// numbers give the Grid values at those coordinates. \p K different Grid%s are
/// returned.
template<typename T, std::size_t N, std::size_t K = 1>
std::array<Grid<T, N>, K> read_grids(
	std::vector<std::array<T, N + K> > const& raw_data,
	T tolerance = 1.e2 * std::numeric_limits<T>::epsilon());

/// Not enough points were provided to form a Grid without ragged edges.
struct NotEnoughPointsError : public std::runtime_error {
	std::size_t points;
	std::size_t expected_points;
	NotEnoughPointsError(std::size_t points, std::size_t expected_points);
};

/// One of the dimensions of a Grid is only one data point thick.
struct SingularDimensionError : public std::runtime_error {
	std::size_t dim;
	SingularDimensionError(std::size_t dim);
};

/// The hyper-cube bounds on a Grid invalid, most likely because the hyper-cube
/// has negative volume.
struct InvalidBoundsError : public std::runtime_error {
	InvalidBoundsError();
};

/// Data points used to construct a Grid are not spaced evenly.
struct InvalidSpacingError : public std::runtime_error {
	InvalidSpacingError();
};

/// Data points used to construct a Grid are provided in an unexpected order.
struct UnexpectedGridPointError : public std::runtime_error {
	std::size_t line_number;
	UnexpectedGridPointError(std::size_t line_number);
};
/// \}

}
}

#include "sidis/extra/interpolate.ipp"

#endif

