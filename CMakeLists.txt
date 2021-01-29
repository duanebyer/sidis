cmake_minimum_required(VERSION 3.1...3.15)
if(${CMAKE_VERSION} VERSION_LESS 3.12)
	cmake_policy(VERSION ${CMAKE_MAJOR_VERSION}.${CMAKE_MINOR_VERSION})
endif()

project(
	Sidis
	VERSION 0.2.0.0
	DESCRIPTION "Event generator for polarized SIDIS with radiative corrections"
	HOMEPAGE_URL "https://github.com/duanebyer/sidis"
	LANGUAGES CXX)

set(
	Sidis_BUILD_TESTS TRUE
	CACHE BOOL
	"should the tests be built (overrides BUILD_TESTING)")
set(
	Sidis_BUILD_EXAMPLES TRUE
	CACHE BOOL
	"should the examples be built")
set(
	Sidis_BUILD_APPS TRUE
	CACHE BOOL
	"should the apps (generator) be built")
set(
	Sidis_REAL_TYPE "double"
	CACHE STRING
	"what default floating point type should be used for real numbers")
set_property(
	CACHE Sidis_REAL_TYPE
	PROPERTY STRINGS
	"float" "double" "long double")

if(NOT CMAKE_PROJECT_NAME STREQUAL PROJECT_NAME)
	message(FATAL_ERROR "Sidis does not support subproject build")
endif()

set(BUILD_TESTING ${Sidis_BUILD_TESTS})
include(CTest)

configure_file(
	"${Sidis_SOURCE_DIR}/cog/const.py.in"
	"${Sidis_BINARY_DIR}/cog/const.py")
configure_file(
	"${Sidis_SOURCE_DIR}/cog/gen_struct.py"
	"${Sidis_BINARY_DIR}/cog/gen_struct.py"
	COPYONLY)

# Check the configuration.
add_subdirectory(check)

# Add external libraries.
add_subdirectory(external/mstwpdf)
add_subdirectory(external/cubature-cpp)

# Set up the library/app builds.
add_subdirectory(src)
if(${Sidis_BUILD_EXAMPLES})
	add_subdirectory(example)
endif()
if(${Sidis_BUILD_APPS})
	add_subdirectory(app/sidisgen)
endif()

# Set up the test builds.
if(${Sidis_BUILD_TESTS})
	enable_testing()
	add_subdirectory(test)
endif()

