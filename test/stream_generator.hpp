#ifndef SIDIS_TEST_STREAM_GENERATOR_HPP
#define SIDIS_TEST_STREAM_GENERATOR_HPP

#include <catch2/catch.hpp>

#include <istream>
#include <memory>
#include <stdexcept>
#include <string>
#include <utility>

template<typename S, typename T>
class StreamGenerator final : public Catch::Generators::IGenerator<T> {
	S _stream;
	T _value;
	std::size_t _elem_count;

public:
	StreamGenerator(S&& stream, bool skip_header = false) :
			_stream(std::move(stream)),
			_elem_count(0) {
		if (skip_header) {
			std::string header;
			std::getline(_stream, header);
			_elem_count += 1;
		}
		if (!next()) {
			throw std::runtime_error(
				"StreamGenerator couldn't read first item from stream");
		}
	}

	bool next() override {
		_stream >> _value;
		_elem_count += 1;
		if (!_stream.eof() && _stream.fail()) {
			throw std::runtime_error(
				"StreamGenerator was unable to parse element "
				+ std::to_string(_elem_count) + " from stream");
		}
		return (bool) _stream;
	}

	T const& get() const override {
		return _value;
	}
};

template<typename T, typename S>
Catch::Generators::GeneratorWrapper<T> from_stream(
		S&& stream,
		bool skip_header = false) {
	return Catch::Generators::GeneratorWrapper<T>(
		std::unique_ptr<Catch::Generators::IGenerator<T> >(
			new StreamGenerator<S, T>(std::move(stream), skip_header)));
}

#endif

