cmake_minimum_required(VERSION 3.1...3.15)
if(${CMAKE_VERSION} VERSION_LESS 3.12)
	cmake_policy(VERSION ${CMAKE_MAJOR_VERSION}.${CMAKE_MINOR_VERSION})
endif()

project(
	Sidis
	VERSION 0.3.0.0
	DESCRIPTION "Event generator for polarized SIDIS with radiative corrections"
	HOMEPAGE_URL "https://github.com/duanebyer/sidis"
	LANGUAGES CXX)
if(NOT CMAKE_PROJECT_NAME STREQUAL PROJECT_NAME)
	message(FATAL_ERROR "Sidis does not support subproject build")
endif()

set(
	Sidis_BUILD_TESTS NO
	CACHE BOOL
	"should the tests be built (overrides BUILD_TESTING)")
set(
	Sidis_BUILD_EXAMPLES YES
	CACHE BOOL
	"should the examples be built")
set(
	Sidis_BUILD_APPS YES
	CACHE BOOL
	"should the apps (generator) be built")

# Documentation is by default built if Latex was found.
find_package(LATEX)
set(
	Sidis_BUILD_DOCS ${LATEX_FOUND}
	CACHE BOOL
	"should the documentation be generated")
# Doxygen documentation is by default built if Doxygen was found.
find_package(Doxygen)
set(
	Sidis_BUILD_DOXYGEN ${DOXYGEN_FOUND}
	CACHE BOOL
	"should the Doxygen documentation be generated")

set(
	Sidis_REAL_TYPE "double"
	CACHE STRING
	"what default floating point type should be used for real numbers")
set_property(
	CACHE Sidis_REAL_TYPE
	PROPERTY STRINGS
	"float" "double" "long double")

set(IPO_SUPPORTED YES)
if(${CMAKE_VERSION} VERSION_GREATER_EQUAL 3.9)
	include(CheckIPOSupported)
	check_ipo_supported(RESULT IPO_SUPPORTED OUTPUT IPO_OUTPUT)
	if(NOT IPO_SUPPORTED)
		message(WARNING "IPO not supported: ${IPO_OUTPUT}")
	endif()
endif()
set(
	Sidis_IPO_ENABLED ${IPO_SUPPORTED}
	CACHE BOOL
	"should interprocess optimization be enabled")

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
add_subdirectory(external/mstwpdf EXCLUDE_FROM_ALL)
add_subdirectory(external/cubature-cpp EXCLUDE_FROM_ALL)

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
# Set up the documentation builds.
if(${Sidis_BUILD_DOCS})
	add_subdirectory(doc)
endif()

