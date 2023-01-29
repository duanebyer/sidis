#ifndef SIDISGEN_EXCEPTION_HPP
#define SIDISGEN_EXCEPTION_HPP

#include <exception>
#include <string>

int const SUCCESS                             = 0;
int const ERROR                               = -1;
int const ERROR_ARG_PARSE                     = -2;
int const ERROR_FILE_NOT_FOUND                = -3;
int const ERROR_FILE_NOT_CREATED              = -4;
int const ERROR_NO_EVENT_TYPES_ENABLED        = -5;
int const ERROR_WRITING_FOAM                  = -6;
int const ERROR_READING_FOAM                  = -7;
int const ERROR_BUILDING_FOAM                 = -8;
int const ERROR_READING_PARAMS                = -9;
int const ERROR_WRITING_PARAMS                = -10;
int const ERROR_PARAMS_PARSE                  = -11;
int const ERROR_PARAMS_INVALID                = -12;
int const ERROR_MERGING_PARAMS                = -13;
int const ERROR_FOAM_INCOMPATIBLE             = -14;
int const ERROR_FOAM_NOT_FOUND                = -15;
int const ERROR_STRUCTURE_FUNCTIONS_NOT_FOUND = -16;
int const ERROR_STRUCTURE_FUNCTIONS_PARSE     = -17;
int const ERROR_WRITING_EVENTS                = -18;
int const ERROR_READING_STATS                 = -19;
int const ERROR_WRITING_STATS                 = -20;
int const ERROR_UNIMPLEMENTED                 = -21;

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

