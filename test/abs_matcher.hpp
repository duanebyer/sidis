#ifndef SIDIS_TEST_ABS_MATCHER_HPP
#define SIDIS_TEST_ABS_MATCHER_HPP

#include <catch2/catch.hpp>

#include <cmath>
#include <iomanip>
#include <ios>
#include <string>
#include <sstream>

template<typename T>
class AbsMatcher final : public Catch::MatcherBase<T> {
	T _value;
	T _abs_error;
	mutable T _cache;

public:
	AbsMatcher(T value, T abs_error) :
		_value(value),
		_abs_error(std::abs(abs_error)) { }

	bool match(T const& x) const override {
		// Using this cache variable is a very bad idea, but it's necessary to
		// work around how `Catch2` deals with printing floating point numbers.
		_cache = x;
		return std::abs(_value - x) <= _abs_error;
	}

	virtual std::string describe() const override {
		int precision =
			(int) std::ceil(std::fmax(-std::log10(_abs_error), 0)) + 2;
		std::ostringstream ss;
		ss << std::scientific << std::setprecision(precision);
		ss << "(" << _cache << ")";
		ss << std::defaultfloat << std::setprecision(3);
		ss << " is within " << _abs_error;
		ss << std::scientific << std::setprecision(precision);
		ss << " of " << _value;
		return ss.str();
	}
};

#endif

