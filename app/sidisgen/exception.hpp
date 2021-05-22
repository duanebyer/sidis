#ifndef SIDISGEN_EXCEPTION_HPP
#define SIDISGEN_EXCEPTION_HPP

#include <exception>
#include <string>

int const SUCCESS = 0;
int const ERROR = -1;
int const ERROR_ARG_PARSE = -2;
int const ERROR_FILE_NOT_FOUND = -3;
int const ERROR_FILE_NOT_CREATED = -4;
int const ERROR_PARAMS_PARSE = -5;
int const ERROR_PARAMS_INVALID = -6;
int const ERROR_FOAM_INCOMPATIBLE = -8;
int const ERROR_FOAM_NOT_FOUND = -9;
int const ERROR_STRUCTURE_FUNCTIONS_NOT_FOUND = -10;
int const ERROR_STRUCTURE_FUNCTIONS_PARSE = -11;
int const ERROR_UNIMPLEMENTED = -12;

class Exception : public std::exception {
	std::string const _what;

public:
	int const error_code;

	Exception(int error_code, char const* what) :
		_what(what),
		error_code(error_code) { }
	Exception(int error_code, std::string what) :
		_what(what),
		error_code(error_code) { }

	char const* what() const noexcept override {
		return _what.c_str();
	}
};

#endif

