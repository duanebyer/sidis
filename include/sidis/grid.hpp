#ifndef SIDIS_GRID_HPP
#define SIDIS_GRID_HPP

#include "sidis/flavor_vec.hpp"
#include "sidis/pimpl.hpp"
#include "sidis/numeric.hpp"

namespace sidis {
namespace sf {
namespace grid {

/**
 * \defgroup GridGroup Loading structure functions, PDFs, and TMDs from grids
 * Utility types used to load data from grids and provide interpolation. These
 * can be used together with SfSet and TmdSet to implement structure functions.
 *
 * Grids can be provided in the LHAPDF, TMDLib, or CSV formats.
 */
/// \{

/**
 * For use with CsvFormat. Documents which variable is stored in which column
 * for a grid stored in a CSV file.
 */
enum class CsvCol {
	/// Column which should be ignored.
	SKIP,
	/// Column with kinematic variable.
	/// \{
	X, Z, Q, Q_SQ, PH_T, PH_T_SQ, QT_TO_Q, K_PERP, K_PERP_SQ,
	/// \}
	/// Column with output (either structure function, PDF value, or TMD value).
	OUT,
};

/**
 * Describes the format of a structure function, PDF, or TMD grid stored in a
 * CSV file. This structure should be filled out for use with CsvSf, CsvPdf,
 * or CsvTmd.
 */
struct CsvFormat {
	/// Number of header rows that should be skipped.
	unsigned header_count = 1;
	/// Number of columns in the CSV file.
	unsigned col_count = 0;
	/// Whether to use a separator between fields. If set to true, the
	/// CsvFormat.separator character is expected between each field. If set to
	/// false, whitespace is used to separate the fields instead.
	bool use_separator = true;
	/// Separator between fields. Ignored if CsvFormat.use_separator is false.
	char separator = ',';
	/// Must be set to an array of length CsvFormat.col_count. Records which
	/// variable is expected to be in each column.
	CsvCol* cols = nullptr;
};

/**
 * Load structure functions from CSV file.
 * \sa CsvFormat
 */
class CsvSf final {
	struct Impl;
	pimpl::Pimpl<Impl> _impl;
public:
	/// Loads structure functions from \p file_name, which must have the format
	/// described by \p format.
	CsvSf(CsvFormat format, char const* file_name);
	/// Evaluate the structure function.
	Real operator()(Real x, Real z, Real Q_sq, Real ph_t_sq) const;
};

/**
 * Load PDFs from CSV file.
 * \sa CsvFormat
 */
struct CsvPdf {
	struct Impl;
	pimpl::Pimpl<Impl> _impl;
public:
	/// Loads PDFs from \p file_name, which must have the format described by
	/// \p format.
	CsvPdf(CsvFormat format, char const* file_name);
	/// Evaluate the PDFs. The flavors are returned in the same order as given
	/// in the CSV.
	FlavorVec operator()(Real x, Real Q_sq) const;
};

/**
 * Load TMDs from CSV file.
 * \sa CsvFormat
 */
struct CsvTmd {
	struct Impl;
	pimpl::Pimpl<Impl> _impl;
public:
	/// Loads TMDs from \p file_name, which must have the format described by
	/// \p format.
	CsvTmd(CsvFormat format, char const* file_name);
	/// Evaluate the TMDs. The flavors are returned in the same order as given
	/// in the CSV.
	FlavorVec operator()(Real x, Real Q_sq, Real k_perp_sq) const;
};

/**
 * Loads and provides access to a PDF from an LHAPDF set.
 */
class Lhapdf final {
	struct Impl;
	pimpl::Pimpl<Impl> _impl;
public:
	/// Loads member \p member from set named \p name.
	Lhapdf(char const* name, int member);
	/// Calculates \f$ x f(x, Q^2) \f$. The resulting FlavorVec is given in the
	/// LHAPDF5 order:
	/// \f$ [\bar{t},\bar{b},\bar{c},\bar{s},\bar{u},\bar{d},g,d,u,s,c,b,t] \f$.
	FlavorVec operator()(Real x, Real Q_sq) const;
	/// Returns quark charges.
	FlavorVec charge() const;
};

struct Tmdlib {
	struct Impl;
	pimpl::Pimpl<Impl> _impl;
public:
	Tmdlib(char const* name, int member);
	FlavorVec operator()(Real x, Real Q_sq, Real k_perp_sq) const;
	FlavorVec charge() const;
};
/// \}

}
}
}

#endif

