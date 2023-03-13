#ifndef SIDIS_FLAVOR_VEC_HPP
#define SIDIS_FLAVOR_VEC_HPP

#include <initializer_list>

#include "sidis/numeric.hpp"

namespace sidis {
namespace sf {

/// Maximum number of TMD flavors supported in user-defined TMDs.
unsigned const MAX_FLAVOR_VEC_SIZE = 15;

/**
 * An aggregate of results from a TMD or FF calculation (from TmdSet), one entry
 * for each parton flavor. The contributions from each flavor are eventually
 * summed together to make a complete structure function by SfSet.
 *
 * Internally, a FlavorVec is a stack-allocated array. Its size cannot be larger
 * than the constant MAX_FLAVOR_ARRAY_SIZE.
 *
 * For convenience, FlavorVec%s can be added, subtracted, and multiplied
 * together.
 *
 * \sa TmdSet
 */
class FlavorVec final {
	// Implemented basically as a smallvec.
	Real _arr[MAX_FLAVOR_VEC_SIZE];
	unsigned _size;
	// Raises an error when a FlavorVec is constructed with a size greater than
	// MAX_FLAVOR_ARRAY_SIZE.
	void raise_out_of_range(unsigned size) const;

public:
	/// Construct a zero-valued FlavorVec of length \p count, which must be less
	/// than MAX_FLAVOR_ARRAY_SIZE.
	explicit FlavorVec(unsigned count) : _arr{}, _size(count) {
		if (_size > MAX_FLAVOR_VEC_SIZE) {
			raise_out_of_range(_size);
		}
	}
	/// Construct a constant-valued FlavorVec of \p count entries of \p value.
	/// \p count must be less than MAX_FLAVOR_ARRAY_SIZE.
	FlavorVec(unsigned count, Real value) : _arr{}, _size(count) {
		if (_size > MAX_FLAVOR_VEC_SIZE) {
			raise_out_of_range(_size);
		}
		for (unsigned fl = 0; fl < _size; ++fl) {
			_arr[fl] = value;
		}
	}
	/// Initialize a FlavorVec from a range.
	template<typename It>
	FlavorVec(It begin, It end) {
		for (It it = begin; it != end; ++it) {
			if (_size >= MAX_FLAVOR_VEC_SIZE) {
				raise_out_of_range(_size + 1);
			}
			_arr[_size] = *it;
			++_size;
		}
	}
	/// Initialize a FlavorVec from an initializer list.
	FlavorVec(std::initializer_list<Real> list);

	/// \name Arithmetic operations
	/// Note that these operations are undefined if the operands have different
	/// sizes (i.e. adding together a length 3 and 4 FlavorVec).
	/// \{
	FlavorVec& operator*=(Real scale) {
		for (unsigned fl = 0; fl < _size; ++fl) {
			_arr[fl] *= scale;
		}
		return *this;
	}
	FlavorVec& operator/=(Real scale) {
		for (unsigned fl = 0; fl < _size; ++fl) {
			_arr[fl] /= scale;
		}
		return *this;
	}
	FlavorVec& operator+=(Real offset) {
		for (unsigned fl = 0; fl < _size; ++fl) {
			_arr[fl] += offset;
		}
		return *this;
	}
	FlavorVec& operator-=(Real offset) {
		for (unsigned fl = 0; fl < _size; ++fl) {
			_arr[fl] -= offset;
		}
		return *this;
	}
	FlavorVec& operator+=(FlavorVec const& rhs) {
		for (unsigned fl = 0; fl < _size; ++fl) {
			_arr[fl] += rhs._arr[fl];
		}
		return *this;
	}
	FlavorVec& operator-=(FlavorVec const& rhs) {
		for (unsigned fl = 0; fl < _size; ++fl) {
			_arr[fl] -= rhs._arr[fl];
		}
		return *this;
	}
	FlavorVec& operator*=(FlavorVec const& rhs) {
		for (unsigned fl = 0; fl < _size; ++fl) {
			_arr[fl] *= rhs._arr[fl];
		}
		return *this;
	}
	FlavorVec& operator/=(FlavorVec const& rhs) {
		for (unsigned fl = 0; fl < _size; ++fl) {
			_arr[fl] /= rhs._arr[fl];
		}
		return *this;
	}
	FlavorVec& inv(Real scale=1.) {
		for (unsigned fl = 0.; fl < _size; ++fl) {
			_arr[fl] = scale / _arr[fl];
		}
		return *this;
	}
	/// \}

	/// Totals together all of the entries.
	Real sum() const {
		Real result = 0.;
		for (unsigned fl = 0; fl < _size; ++fl) {
			result += _arr[fl];
		}
		return result;
	}

	/// \name Array access
	/// \{
	unsigned size() const {
		return _size;
	}
	Real const& operator[](unsigned fl) const {
		return _arr[fl];
	}
	Real& operator[](unsigned fl) {
		return _arr[fl];
	}
	Real* data() {
		return _arr;
	}
	Real const* data() const {
		return _arr;
	}
	/// \}

	friend FlavorVec tmd_gaussian_factor(FlavorVec var, Real k_perp_sq);
	friend inline FlavorVec operator-(FlavorVec vec) {
		for (unsigned fl = 0; fl < vec._size; ++fl) {
			vec._arr[fl] = -vec._arr[fl];
		}
		return vec;
	}
};

// Arithmetic with real numbers.
inline FlavorVec operator*(FlavorVec lhs, Real rhs) {
	lhs *= rhs;
	return lhs;
}
inline FlavorVec operator*(Real lhs, FlavorVec rhs) {
	rhs *= lhs;
	return rhs;
}

inline FlavorVec operator/(FlavorVec lhs, Real rhs) {
	lhs /= rhs;
	return lhs;
}
inline FlavorVec operator/(Real lhs, FlavorVec rhs) {
	rhs.inv(lhs);
	return rhs;
}

inline FlavorVec operator+(FlavorVec lhs, Real rhs) {
	lhs += rhs;
	return lhs;
}
inline FlavorVec operator+(Real lhs, FlavorVec rhs) {
	rhs += lhs;
	return rhs;
}

inline FlavorVec operator-(FlavorVec lhs, Real rhs) {
	lhs -= rhs;
	return lhs;
}
inline FlavorVec operator-(Real lhs, FlavorVec rhs) {
	rhs -= lhs;
	return -rhs;
}

// Arithemetic with two FlavorVec%s.
inline FlavorVec operator+(FlavorVec lhs, FlavorVec const& rhs) {
	lhs += rhs;
	return lhs;
}
inline FlavorVec operator-(FlavorVec lhs, FlavorVec const& rhs) {
	lhs -= rhs;
	return lhs;
}
inline FlavorVec operator*(FlavorVec lhs, FlavorVec const& rhs) {
	lhs *= rhs;
	return lhs;
}
inline FlavorVec operator/(FlavorVec lhs, FlavorVec const& rhs) {
	lhs /= rhs;
	return lhs;
}

FlavorVec sqrt_vec(FlavorVec vec);
FlavorVec sq_vec(FlavorVec vec);
FlavorVec exp_vec(FlavorVec vec);
FlavorVec pow_vec(Real lhs, FlavorVec rhs);
FlavorVec pow_vec(FlavorVec lhs, Real rhs);
FlavorVec pow_vec(FlavorVec lhs, FlavorVec const& rhs);

/// Convenience function for computing the factor that shows up in many
/// Gaussian TMDs. The flavor-dependent variance of the Gaussian is given by
/// \p var, and the momentum at which to evaluate the Gaussian is \p k_perp_sq.
FlavorVec tmd_gaussian_factor(FlavorVec var, Real k_perp_sq);

}
}

#endif

