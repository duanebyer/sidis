#include "sidis/grid.hpp"

#include <algorithm>
#include <array>
#include <fstream>
#include <limits>
#include <memory>
#include <utility>
#include <type_traits>
#include <vector>

#include <LHAPDF/LHAPDF.h>

#include "sidis/utility.hpp"
#include "sidis/extra/exception.hpp"
#include "sidis/extra/interpolate.hpp"

using namespace sidis;
using namespace sidis::interp;
using namespace sidis::sf;
using namespace sidis::sf::grid;

#define SF_SET_DIR "sidis/sf_set"

namespace {

FlavorVec const LHAPDF_CHARGE = {
	+1./3., -2./3., +1./3., -2./3., +1./3., -2./3., // anti-quarks: t b c s u d
	0.,                                             // gluons
	-1./3., +2./3., -1./3., +2./3., -1./3., +2./3., // quarks:      d u s c b t
};

// Using this trick to statically initialize the path prefixes one time at the
// beginning of the program.
static bool const set_lhapdf_path_prefixes = []() {
	LHAPDF::pathsPrepend(DATADIR "/" SF_SET_DIR "/");
	LHAPDF::pathsPrepend("../share/" SF_SET_DIR "/");
	LHAPDF::pathsPrepend(SF_SET_DIR "/");
	LHAPDF::pathsPrepend("./");
	return true;
}();

// Searches several standard directories for a file, and opens it once found. If
// the file can't be found, throws an exception.
std::ifstream& find_file(std::ifstream& fin, char const* file_name) {
	// TODO: This path lookup code is very broken, it doesn't work correctly
	// when the provided file name is an absolute path, or on operating systems
	// with a different path separator. It would be good to replace it with a
	// path handling library in the future.
	fin.open(std::string(DATADIR "/" SF_SET_DIR "/") + file_name);
	if (fin) {
		return fin;
	}
	fin.open(std::string("../share/" SF_SET_DIR "/") + file_name);
	if (fin) {
		return fin;
	}
	fin.open(std::string(SF_SET_DIR "/") + file_name);
	if (fin) {
		return fin;
	}
	fin.open(file_name);
	if (fin) {
		return fin;
	}
	fin.open(file_name);
	if (fin) {
		return fin;
	}
	throw DataFileNotFound(file_name);
}

std::size_t INVALID_INDEX = std::numeric_limits<std::size_t>::max();

enum class CsvType {
	SF,
	PDF,
	TMD,
};

// Implements the CSV grid loading and lookup.
template<CsvType CT>
struct CsvImpl {
	static constexpr std::size_t SLOT_COUNT =
		CT == CsvType::SF ? 4
		: CT == CsvType::PDF ? 2
		: CT == CsvType::TMD ? 3
		: 0;
	static std::size_t axis_assign(CsvCol col);

	// The grids themselves, all with the same dimensions, sizes, and bounds.
	std::vector<Grid<Real, SLOT_COUNT> > grids;
	// The interpolated view into the grids.
	std::vector<CubicView<Real, SLOT_COUNT> > views;
	// The transformation for each axis.
	std::array<CsvCol, SLOT_COUNT> axes_types;

	CsvImpl(CsvFormat format, char const* file_name);
	Real operator()(std::size_t fl, std::array<Real, SLOT_COUNT> x) const;
};

template<CsvType CT>
CsvImpl<CT>::CsvImpl(CsvFormat format, char const* file_name) {
	// Find a map from CSV indices (called "columns") to grid indices (called
	// "axes"). The map goes both ways.
	std::vector<std::size_t> col_to_axis(format.col_count, INVALID_INDEX);
	std::vector<std::size_t> axis_to_col(SLOT_COUNT, INVALID_INDEX);
	// Similarly, find a bidirectional map from CSV indices ("columns") to
	// flavor indices.
	std::vector<std::size_t> col_to_fl(format.col_count, INVALID_INDEX);
	std::vector<std::size_t> fl_to_col;
	if (format.cols == nullptr) {
		throw CsvFormatMalformed(file_name);
	}
	for (std::size_t col_idx = 0; col_idx < format.col_count; ++col_idx) {
		CsvCol col = format.cols[col_idx];
		if (col == CsvCol::OUT) {
			// It's an output column. Add it to the list.
			col_to_fl[col_idx] = fl_to_col.size();
			fl_to_col.push_back(col_idx);
		} else {
			// It's an input column. Store it in the map.
			std::size_t axis_idx = axis_assign(col);
			col_to_axis[col_idx] = axis_idx;
			if (axis_idx != INVALID_INDEX) {
				if (axis_to_col[axis_idx] != INVALID_INDEX) {
					// Error: this axis has already been assigned a column.
					throw CsvFormatMalformed(file_name);
				}
				axis_to_col[axis_idx] = col_idx;
			}
		}
	}
	for (std::size_t axis_idx = 0; axis_idx < SLOT_COUNT; ++axis_idx) {
		if (axis_to_col[axis_idx] == INVALID_INDEX) {
			// Error: this axis didn't get assigned to any column.
			throw CsvFormatMalformed(file_name);
		}
		axes_types[axis_idx] = format.cols[axis_to_col[axis_idx]];
	}
	std::size_t fl_count = fl_to_col.size();
	std::size_t stride = SLOT_COUNT + fl_count;
	std::size_t row_count = 0;

	// Read the file into a data vector.
	std::ifstream in;
	find_file(in, file_name);
	std::vector<Real> data;
	if (format.use_separator) {
		std::string line;
		while (std::getline(in, line)) {
			std::stringstream line_in(line);
			// Add the next column to the data vector.
			data.resize(
				data.size() + stride,
				std::numeric_limits<Real>::quiet_NaN());
			auto row = data.begin() + stride * row_count;
			for (std::size_t col_idx = 0; col_idx < format.col_count; ++col_idx) {
				if (format.use_separator && col_idx > 0) {
					char separator;
					line_in >> separator;
					if (separator != format.separator) {
						throw DataFileParseError(file_name);
					}
				}
				Real field;
				line_in >> field;
				if (col_to_fl[col_idx] != INVALID_INDEX) {
					row[SLOT_COUNT + col_to_fl[col_idx]] = field;
				} else if (col_to_axis[col_idx] != INVALID_INDEX) {
					row[col_to_axis[col_idx]] = field;
				}
			}
			row_count += 1;
			if (!line_in) {
				throw DataFileParseError(file_name);
			}
		}
	}

	// By default, assume the grids are only accurate to single precision.
	grids = read_grids<Real, SLOT_COUNT>(data, fl_count, 0.000001);

	// Set up the views to look into the grids.
	views.reserve(fl_count);
	for (std::size_t fl_idx = 0; fl_idx < fl_count; ++fl_idx) {
		views.emplace_back(grids[fl_idx]);
	}
}

// These `axis_assign` functions take the type of variable in a column and
// assign it to specific index (e.g. Q -> 2). This is used to ensure that if the
// columns of a CSV file are given in a different order than expected, the
// variables can be re-ordered into a standard order.
template<>
std::size_t CsvImpl<CsvType::SF>::axis_assign(CsvCol col) {
	switch (col) {
	case CsvCol::X:
		return 0;
	case CsvCol::Z:
		return 1;
	case CsvCol::Q:
	case CsvCol::Q_SQ:
		return 2;
	case CsvCol::PH_T:
	case CsvCol::PH_T_SQ:
	case CsvCol::QT_TO_Q:
		return 3;
	default:
		return INVALID_INDEX;
	}
}
template<>
std::size_t CsvImpl<CsvType::PDF>::axis_assign(CsvCol col) {
	switch (col) {
	case CsvCol::X:
		return 0;
	case CsvCol::Q:
	case CsvCol::Q_SQ:
		return 1;
	default:
		return INVALID_INDEX;
	}
}
template<>
std::size_t CsvImpl<CsvType::TMD>::axis_assign(CsvCol col) {
	switch (col) {
	case CsvCol::X:
		return 0;
	case CsvCol::Q:
	case CsvCol::Q_SQ:
		return 1;
	case CsvCol::K_PERP:
	case CsvCol::K_PERP_SQ:
		return 2;
	default:
		return INVALID_INDEX;
	}
}

// Evaluate the grids with cubic interpolation. The input kinematics must be
// provided in a standard form (e.g. for structure functions, `x`, `z`, `Q_sq`,
// and `ph_t_sq`). The kinematics are then transformed into the internal
// coordinates used by the grid. For example, the grid may use `Q` instead of
// `Q_sq`.
template<>
Real CsvImpl<CsvType::SF>::operator()(std::size_t fl, std::array<Real, 4> p) const {
	Real x = p[0];
	Real z = p[1];
	Real Q_sq = p[2];
	Real ph_t_sq = p[3];
	switch (axes_types[0]) {
	case CsvCol::X:
		p[0] = x;
		break;
	default:
		SIDIS_UNREACHABLE();
	}
	switch (axes_types[1]) {
	case CsvCol::Z:
		p[1] = z;
		break;
	default:
		SIDIS_UNREACHABLE();
	}
	switch (axes_types[2]) {
	case CsvCol::Q:
		p[2] = std::sqrt(Q_sq);
		break;
	case CsvCol::Q_SQ:
		p[2] = Q_sq;
		break;
	default:
		SIDIS_UNREACHABLE();
	}
	switch (axes_types[3]) {
	case CsvCol::PH_T:
		p[3] = std::sqrt(ph_t_sq);
		break;
	case CsvCol::PH_T_SQ:
		p[3] = ph_t_sq;
		break;
	case CsvCol::QT_TO_Q:
		p[3] = std::sqrt(ph_t_sq / Q_sq) / z;
		break;
	default:
		SIDIS_UNREACHABLE();
	}
	return views[fl](p);
}
template<>
Real CsvImpl<CsvType::PDF>::operator()(std::size_t fl, std::array<Real, 2> p) const {
	Real x = p[0];
	Real Q_sq = p[1];
	switch (axes_types[0]) {
	case CsvCol::X:
		p[0] = x;
		break;
	default:
		SIDIS_UNREACHABLE();
	}
	switch (axes_types[1]) {
	case CsvCol::Q:
		p[1] = std::sqrt(Q_sq);
		break;
	case CsvCol::Q_SQ:
		p[1] = Q_sq;
		break;
	default:
		SIDIS_UNREACHABLE();
	}
	return views[fl](p);
}
template<>
Real CsvImpl<CsvType::TMD>::operator()(std::size_t fl, std::array<Real, 3> p) const {
	Real x = p[0];
	Real Q_sq = p[1];
	Real k_perp_sq = p[2];
	switch (axes_types[0]) {
	case CsvCol::X:
		p[0] = x;
		break;
	default:
		SIDIS_UNREACHABLE();
	}
	switch (axes_types[1]) {
	case CsvCol::Q:
		p[1] = std::sqrt(Q_sq);
		break;
	case CsvCol::Q_SQ:
		p[1] = Q_sq;
		break;
	default:
		SIDIS_UNREACHABLE();
	}
	switch (axes_types[2]) {
	case CsvCol::K_PERP:
		p[2] = std::sqrt(k_perp_sq);
		break;
	case CsvCol::K_PERP_SQ:
		p[2] = k_perp_sq;
		break;
	default:
		SIDIS_UNREACHABLE();
	}
	return views[fl](p);
}

}

struct CsvSf::Impl : public CsvImpl<CsvType::SF> {
	Impl(CsvFormat format, char const* file_name) : CsvImpl(format, file_name) { }
};
struct CsvPdf::Impl : public CsvImpl<CsvType::PDF> {
	Impl(CsvFormat format, char const* file_name) : CsvImpl(format, file_name) { }
};
struct CsvTmd::Impl : public CsvImpl<CsvType::TMD> {
	Impl(CsvFormat format, char const* file_name) : CsvImpl(format, file_name) { }
};
struct Lhapdf::Impl {
	std::unique_ptr<LHAPDF::PDF> pdf_set;
	Impl(char const* name, int member) :
		pdf_set(LHAPDF::mkPDF(name, member)) { }
};
struct Tmdlib::Impl {
	Impl() { }
};

CsvSf::CsvSf(CsvFormat format, char const* file_name) :
	_impl(pimpl::make_pimpl<Impl>(format, file_name)) { }
CsvPdf::CsvPdf(CsvFormat format, char const* file_name) :
	_impl(pimpl::make_pimpl<Impl>(format, file_name)) { }
CsvTmd::CsvTmd(CsvFormat format, char const* file_name) :
	_impl(pimpl::make_pimpl<Impl>(format, file_name)) { }
Lhapdf::Lhapdf(char const* name, int member) :
	_impl(pimpl::make_pimpl<Impl>(name, member)) { }

Real CsvSf::operator()(Real x, Real z, Real Q_sq, Real ph_t_sq) const {
	return (*_impl)(0, { x, z, Q_sq, ph_t_sq });
}
FlavorVec CsvPdf::operator()(Real x, Real Q_sq) const {
	FlavorVec result(_impl->views.size());
	for (unsigned fl = 0; fl < _impl->views.size(); ++fl) {
		result[fl] = (*_impl)(fl, { x, Q_sq });
	}
	return result;
}
FlavorVec CsvTmd::operator()(Real x, Real Q_sq, Real k_perp_sq) const {
	FlavorVec result(_impl->views.size());
	for (unsigned fl = 0; fl < _impl->views.size(); ++fl) {
		result[fl] = (*_impl)(fl, { x, Q_sq, k_perp_sq });
	}
	return result;
}
FlavorVec Lhapdf::operator()(Real x, Real Q_sq) const {
	thread_local std::vector<Real> pdf_out(LHAPDF_CHARGE.size());
	FlavorVec result(LHAPDF_CHARGE.size());
	_impl->pdf_set->xfxQ2(x, Q_sq, pdf_out);
	std::copy(pdf_out.begin(), pdf_out.end(), result.data());
	return result;
}

FlavorVec Lhapdf::charge() const {
	return LHAPDF_CHARGE;
}

