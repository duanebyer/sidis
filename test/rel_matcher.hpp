#ifndef SIDIS_TEST_REL_MATCHER_HPP
#define SIDIS_TEST_REL_MATCHER_HPP

#include <catch2/catch.hpp>

#include <cmath>
#include <iomanip>
#include <ios>
#include <string>
#include <sstream>

template<typename T>
class RelMatcher final : public Catch::MatcherBase<T> {
	T _value;
	T _rel_error;
	mutable T _cache;

public:
	RelMatcher(T value, T rel_error) :
		_value(value),
		_rel_error(std::abs(rel_error)) { }

	bool match(T const& x) const override {
		// Using this cache variable is a very bad idea, but it's necessary to
		// work around how `Catch2` deals with printing floating point numbers.
		_cache = x;
		return std::abs(_value - x) <= std::fmax(
			std::abs(_rel_error * _value),
			std::abs(_rel_error * x));
	}

	virtual std::string describe() const override {
		int precision =
			(int) std::ceil(std::fmax(-std::log10(_rel_error), 0)) + 2;
		std::ostringstream ss;
		ss << std::scientific << std::setprecision(precision);
		ss << "(" << _cache << ")";
		ss << std::defaultfloat << std::setprecision(3);
		ss << " is within fraction " << _rel_error;
		ss << std::scientific << std::setprecision(precision);
		ss << " of " << _value;
		return ss.str();
	}
};

#endif

