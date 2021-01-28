#ifndef SIDIS_EXCEPTION_HPP
#define SIDIS_EXCEPTION_HPP

#include <exception>
#include <string>

#include "sidis/constant.hpp"
#include "sidis/numeric.hpp"

namespace sidis {

/**
 * Mass threshold does not satisfy `M <= Mth <= M + mh`.
 */
class MassThresholdOutOfRange final : public std::exception {
	std::string _what;

public:
	Real Mth;
	Real Mth_min;
	Real Mth_max;

	MassThresholdOutOfRange(Real Mth, Real Mth_min, Real Mth_max);
	char const* what() const noexcept override {
		return _what.c_str();
	}
};

/**
 * Center-of-mass energy is below threshold.
 */
class ComEnergyOutOfRange final : public std::exception {
	std::string _what;

public:
	Real E_b;
	Real E_b_min;

	ComEnergyOutOfRange(Real E_b, Real E_b_min);
	char const* what() const noexcept override {
		return _what.c_str();
	}
};

/**
 * Target nucleus is a different type than expected.
 */
class TargetMismatch final : public std::exception {
	std::string _what;

public:
	constant::Nucleus target;
	constant::Nucleus target_expected;

	TargetMismatch(constant::Nucleus target, constant::Nucleus target_expected);
	char const* what() const noexcept override {
		return _what.c_str();
	}
};

/**
 * Disallowed flavor index in a TMD or FF.
 *
 * This exception is intended for use in user-defined TMDs and FFs to indicate
 * an out-of-bound flavor.
 */
class FlavorOutOfRange final : public std::exception {
	std::string _what;

public:
	unsigned flavor;

	FlavorOutOfRange(unsigned flavor);
	char const* what() const noexcept override {
		return _what.c_str();
	}
};

/**
 * Disallowed hadron in a FF.
 *
 * This exception is intended for use in user-defined FFs to indicate an
 * unsupported type of hadron.
 */
class HadronOutOfRange final : public std::exception {
	std::string _what;

public:
	constant::Hadron hadron;

	HadronOutOfRange(constant::Hadron hadron);
	char const* what() const noexcept override {
		return _what.c_str();
	}
};

/**
 * Data file not able to be opened.
 */
class DataFileNotFound final : public std::exception {
	std::string _what;

public:
	std::string file_name;

	DataFileNotFound(std::string file_name);
	char const* what() const noexcept override {
		return _what.c_str();
	}
};

/**
 * Data could not be parsed.
 */
class DataFileParseError final : public std::exception {
	std::string _what;

public:
	std::string file_name;

	DataFileParseError(std::string file_name);
	char const* what() const noexcept override {
		return _what.c_str();
	}
};

}

#endif

