cmake_minimum_required(VERSION 3.1...3.15)
if(${CMAKE_VERSION} VERSION_LESS 3.12)
	cmake_policy(VERSION ${CMAKE_MAJOR_VERSION}.${CMAKE_MINOR_VERSION})
endif()

project(
	Sidis
	VERSION 0.0.1
	DESCRIPTION "Event generator for polarized SIDIS with radiative corrections"
	HOMEPAGE_URL "https://github.com/duanebyer/sidis"
	LANGUAGES CXX)

if(NOT CMAKE_PROJECT_NAME STREQUAL PROJECT_NAME)
	message(FATAL_ERROR "SIDIS does not support subproject build")
endif()

set(
	CMAKE_MODULE_PATH
	${CMAKE_MODULE_PATH}
	"${Sidis_SOURCE_DIR}/cmake/Mathematica")

add_subdirectory(src)
add_subdirectory(app/sidisgen)

enable_testing()
add_subdirectory(test)

