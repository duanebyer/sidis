#ifndef SIDIS_TEST_REL_MATCHER_HPP
#define SIDIS_TEST_REL_MATCHER_HPP

#include <catch2/catch.hpp>

#include <cmath>
#include <iomanip>
#include <ios>
#include <string>
#include <sstream>

template<typename T>
struct RelMatcher : public Catch::MatcherBase<T> {
	T value;
	T rel_error;
	mutable T cache;

	RelMatcher(T value, T rel_error) : value(value), rel_error(rel_error) { }

	bool match(T const& x) const override {
		// Using this cache variable is a very bad idea, but it's necessary to
		// work around how `Catch2` deals with printing floating point numbers.
		cache = x;
		return std::abs(value - x) <= std::abs(rel_error * value);
	}

	virtual std::string describe() const override {
		int precision = (int) std::ceil(std::fmax(-std::log10(rel_error), 0)) + 2;
		std::ostringstream ss;
		ss << std::scientific << std::setprecision(precision);
		ss << "(" << cache << ")";
		ss << std::defaultfloat << std::setprecision(3);
		ss << " is within " << rel_error;
		ss << std::scientific << std::setprecision(precision);
		ss << " of " << value;
		return ss.str();
	}
};

#endif

