configure_file(
	"${Sidis_SOURCE_DIR}/check/check_numeric.cpp.in"
	"${Sidis_BINARY_DIR}/check/check_numeric.cpp")
try_compile(
	CHECK_NUMERIC
	"${Sidis_BINARY_DIR}/check"
	SOURCES "${Sidis_BINARY_DIR}/check/check_numeric.cpp"
	CXX_STANDARD 11
	CXX_EXTENSIONS OFF)
if(NOT CHECK_NUMERIC)
	message(SEND_ERROR "Sidis real type `${Sidis_REAL_TYPE}` is not a supported floating point type")
endif()

if(Sidis_BUILD_TESTS)
	try_compile(
		CHECK_TEST_FLOAT
		"${Sidis_BINARY_DIR}/check"
		SOURCES "${Sidis_SOURCE_DIR}/check/check_test_float.cpp"
		CXX_STANDARD 11
		CXX_EXTENSIONS OFF)
	try_compile(
		CHECK_TEST_DOUBLE
		"${Sidis_BINARY_DIR}/check"
		SOURCES "${Sidis_SOURCE_DIR}/check/check_test_double.cpp"
		CXX_STANDARD 11
		CXX_EXTENSIONS OFF)
	try_compile(
		CHECK_TEST_LONG_DOUBLE
		"${Sidis_BINARY_DIR}/check"
		SOURCES "${Sidis_SOURCE_DIR}/check/check_test_long_double.cpp"
		CXX_STANDARD 11
		CXX_EXTENSIONS OFF)
	configure_file(
		"${Sidis_SOURCE_DIR}/check/check_test_real.cpp.in"
		"${Sidis_BINARY_DIR}/check/check_test_real.cpp")
	try_compile(
		CHECK_TEST_REAL
		"${Sidis_BINARY_DIR}/check"
		SOURCES "${Sidis_BINARY_DIR}/check/check_test_real.cpp"
		CXX_STANDARD 11
		CXX_EXTENSIONS OFF)
	if(NOT CHECK_TEST_FLOAT)
		message(WARNING "Sidis tests do not support compiler `float` type")
	endif()
	if(NOT CHECK_TEST_DOUBLE)
		message(WARNING "Sidis tests do not support compiler `double` type")
	endif()
	if(NOT CHECK_TEST_LONG_DOUBLE)
		message(WARNING "Sidis tests do not support compiler `long double` type")
	endif()
	if(NOT CHECK_TEST_REAL)
		message(WARNING "Sidis tests do not support Sidis real type of `${Sidis_REAL_TYPE}` (an IEEE 754 64-bit interchange type is the minimum requirement for testing)")
	endif()
endif()

