#ifndef SIDIS_SF_MODEL_WW_HPP
#define SIDIS_SF_MODEL_WW_HPP

#include <stdexcept>
#include <string>

#include "sidis/numeric.hpp"
#include "sidis/structure_function.hpp"

namespace sidis {
namespace sf {
namespace model {

/**
 * Wandzura-Wilczek approximation to the structure functions [2].
 */
class WW final {
	struct Impl;
	Impl* _impl;

public:
	WW();
	~WW();

	WW(WW const& other) = delete;
	WW(WW&& other) noexcept;
	WW& operator=(WW const& other) = delete;
	WW& operator=(WW&& other) noexcept;

	SfUU sf_uu(Real x, Real z, Real Q_sq, Real ph_t) const;
	SfUL sf_ul(Real x, Real z, Real Q_sq, Real ph_t) const;
	SfUT sf_ut(Real x, Real z, Real Q_sq, Real ph_t) const;
	SfLU sf_lu(Real x, Real z, Real Q_sq, Real ph_t) const;
	SfLL sf_ll(Real x, Real z, Real Q_sq, Real ph_t) const;
	SfLT sf_lt(Real x, Real z, Real Q_sq, Real ph_t) const;

	SfXU sf_xu(Real x, Real z, Real Q_sq, Real ph_t) const {
		return { sf_uu(x, z, Q_sq, ph_t), sf_lu(x, z, Q_sq, ph_t) };
	}
	SfXL sf_xl(Real x, Real z, Real Q_sq, Real ph_t) const {
		return { sf_ul(x, z, Q_sq, ph_t), sf_ll(x, z, Q_sq, ph_t) };
	}
	SfXT sf_xt(Real x, Real z, Real Q_sq, Real ph_t) const {
		return { sf_ut(x, z, Q_sq, ph_t), sf_lt(x, z, Q_sq, ph_t) };
	}
	SfUX sf_ux(Real x, Real z, Real Q_sq, Real ph_t) const {
		return {
			sf_uu(x, z, Q_sq, ph_t),
			sf_ul(x, z, Q_sq, ph_t),
			sf_ut(x, z, Q_sq, ph_t),
		};
	}
	SfLX sf_lx(Real x, Real z, Real Q_sq, Real ph_t) const {
		return {
			sf_lu(x, z, Q_sq, ph_t),
			sf_ll(x, z, Q_sq, ph_t),
			sf_lt(x, z, Q_sq, ph_t),
		};
	}

	Sf sf(Real x, Real z, Real Q_sq, Real ph_t) const {
		return { sf_ux(x, z, Q_sq, ph_t), sf_lx(x, z, Q_sq, ph_t) };
	}
};

struct DataFileNotFoundException : std::runtime_error {
	std::string file_name;
	DataFileNotFoundException(std::string file_name) :
		std::runtime_error("Couldn't find data file " + file_name),
		file_name(file_name) { }
};

struct DataFileFormatException : std::runtime_error {
	std::string file_name;
	DataFileFormatException(std::string file_name) :
		std::runtime_error("Invalid format for data file " + file_name),
		file_name(file_name) { }
};

}
}
}

#endif

