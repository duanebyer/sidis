cmake_minimum_required(VERSION 3.2...3.15)
if(CMAKE_VERSION VERSION_LESS 3.12)
	cmake_policy(VERSION ${CMAKE_MAJOR_VERSION}.${CMAKE_MINOR_VERSION})
endif()

project(
	Sidis
	VERSION 0.6.0.0
	DESCRIPTION "Event generator for polarized SIDIS with radiative corrections"
	HOMEPAGE_URL "https://github.com/duanebyer/sidis"
	LANGUAGES CXX)
if(NOT CMAKE_PROJECT_NAME STREQUAL PROJECT_NAME)
	message(FATAL_ERROR "Sidis does not support subproject build")
endif()

list(APPEND CMAKE_MODULE_PATH "${CMAKE_SOURCE_DIR}/cmake")

# Load packages. These aren't used yet, but whether they are found or not
# determines defaults for some options.
find_package(Doxygen)
find_package(LATEX)
find_package(OpenMP 3.0)

# Project-specific cache variable for shared library.
if(DEFINED Sidis_BUILD_SHARED_LIBS)
	set(Sidis_BUILD_SHARED_LIBS ${Sidis_BUILD_SHARED_LIBS}
		CACHE BOOL
		"should the library be built as shared (overrides BUILD_SHARED_LIBS")
	set(BUILD_SHARED_LIBS ${Sidis_BUILD_SHARED_LIBS})
else()
	set(Sidis_BUILD_SHARED_LIBS ${BUILD_SHARED_LIBS})
endif()
# Project-specific cache variable for tests.
if(DEFINED Sidis_BUILD_TESTS)
	set(
		Sidis_BUILD_TESTS ${Sidis_BUILD_TESTS}
		CACHE BOOL
		"should the tests be built (overrides BUILD_TESTING)")
	set(BUILD_TESTING ${Sidis_BUILD_TESTS})
else()
	set(Sidis_BUILD_TESTS ${BUILD_TESTING})
endif()

set(
	Sidis_BUILD_EXAMPLES YES
	CACHE BOOL
	"should the examples be built")
set(
	Sidis_BUILD_APPS YES
	CACHE BOOL
	"should the apps (generator) be built")
set(
	Sidis_BUILD_SCRATCH NO
	CACHE BOOL
	"should the scratch be built")

# Documentation is by default built if Latex was found.
set(
	Sidis_BUILD_DOCS ${LATEX_FOUND}
	CACHE BOOL
	"should the documentation be generated")
# Doxygen documentation is by default built if Doxygen was found.
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

# Interprocedural optimization support.
set(IPO_SUPPORTED NO)
if(CMAKE_VERSION VERSION_GREATER_EQUAL 3.9)
	include(CheckIPOSupported)
	check_ipo_supported(RESULT IPO_SUPPORTED OUTPUT IPO_OUTPUT)
	if(NOT IPO_SUPPORTED)
		message(WARNING "IPO not supported: ${IPO_OUTPUT}")
	endif()
endif()
if(Sidis_BUILD_SHARED_LIBS)
	set(IPO_SUPPORTED NO)
endif()
set(
	Sidis_ENABLE_IPO ${IPO_SUPPORTED}
	CACHE BOOL
	"should interprocedural optimization support be enabled")

# OpenMP multithreading support.
set(OPENMP_SUPPORTED ${OPENMP_FOUND})
if(NOT OpenMP_FOUND)
	message(WARNING "OpenMP not found, multi-threading support disabled")
endif()
set(
	Sidis_OPENMP_ENABLED ${OPENMP_SUPPORTED}
	CACHE BOOL
	"should OpenMP multithreading support be enabled")

# Set up warnings.
set(
	Sidis_COMPILER_WARNINGS
	$<$<OR:$<CXX_COMPILER_ID:AppleClang>,$<CXX_COMPILER_ID:Clang>,$<CXX_COMPILER_ID:GNU>>:-Wall -Wextra -pedantic>
	$<$<CXX_COMPILER_ID:MSVC>:/W4>)

# Add external libraries.
add_subdirectory(external/bubble EXCLUDE_FROM_ALL)
add_subdirectory(external/cubature-cpp EXCLUDE_FROM_ALL)
add_subdirectory(external/mstwpdf EXCLUDE_FROM_ALL)

# Check the configuration.
add_subdirectory(check)

include(CTest)

# Set up the library/app builds.
add_subdirectory(src)
if(Sidis_BUILD_EXAMPLES)
	add_subdirectory(example)
endif()
if(Sidis_BUILD_APPS)
	add_subdirectory(app/sidisgen)
endif()
if(Sidis_BUILD_SCRATCH)
	add_subdirectory(scratch)
endif()
# Set up the test builds.
if(Sidis_BUILD_TESTS)
	enable_testing()
	add_subdirectory(test)
endif()
# Set up the documentation builds.
if(Sidis_BUILD_DOCS)
	add_subdirectory(doc)
endif()

