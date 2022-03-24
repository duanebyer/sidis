set(
	HEADER_LIST_BASE
	"sidis/sidis.hpp"
	"sidis/asymmetry.hpp"
	"sidis/bound.hpp"
	"sidis/constant.hpp"
	"sidis/cross_section.hpp"
	"sidis/cut.hpp"
	"sidis/frame.hpp"
	"sidis/integ.hpp"
	"sidis/kinematics.hpp"
	"sidis/particle.hpp"
	"sidis/tmd.hpp"
	"sidis/transform.hpp"
	"sidis/vector.hpp"
	"sidis/sf_set/mask.hpp"
	"sidis/sf_set/prokudin.hpp"
	"sidis/sf_set/stf.hpp"
	"sidis/sf_set/test.hpp"
	"sidis/extra/exception.hpp"
	"sidis/extra/integrate.hpp"
	"sidis/extra/interpolate.hpp"
	"sidis/extra/interpolate.ipp"
	"sidis/extra/math.hpp"
	"sidis/extra/map.hpp")
set(
	HEADER_LIST_GENERATE
	"sidis/hadronic_coeff.hpp"
	"sidis/leptonic_coeff.hpp"
	"sidis/numeric.hpp"
	"sidis/structure_function.hpp"
	"sidis/version.hpp")

set(HEADER_LIST_SOURCE)
set(HEADER_LIST_DIST)
foreach(HEADER ${HEADER_LIST_BASE})
	list(APPEND HEADER_LIST_SOURCE "${Sidis_SOURCE_DIR}/include/${HEADER}")
	list(APPEND HEADER_LIST_DIST "${Sidis_SOURCE_DIR}/include/${HEADER}")
endforeach(HEADER)
foreach(HEADER ${HEADER_LIST_GENERATE})
	list(APPEND HEADER_LIST_SOURCE "${Sidis_SOURCE_DIR}/include/${HEADER}.in")
	list(APPEND HEADER_LIST_DIST "${Sidis_BINARY_DIR}/include/${HEADER}")
endforeach(HEADER)

# Check that Python is installed.
if(CMAKE_VERSION VERSION_LESS 3.12)
	find_package(PythonInterp 2.7 REQUIRED)
	find_package(PythonLibs 2.7 REQUIRED)
	set(Python_EXECUTABLE PYTHON_EXECUTABLE)
else()
	find_package(Python REQUIRED COMPONENTS Interpreter Development)
endif()

# Generated header files are processed through the `cog` utility.
foreach(HEADER ${HEADER_LIST_GENERATE})
	add_custom_command(
		OUTPUT "${Sidis_BINARY_DIR}/include/${HEADER}"
		COMMAND ${CMAKE_COMMAND}
		ARGS
			-E env "PYTHONPATH=${Sidis_SOURCE_DIR}/external"
			"${Python_EXECUTABLE}" ARGS -m cogapp -d
			-I "${Sidis_BINARY_DIR}/cog"
			-o "${Sidis_BINARY_DIR}/include/${HEADER}"
			"${Sidis_SOURCE_DIR}/include/${HEADER}.in"
		MAIN_DEPENDENCY "${Sidis_SOURCE_DIR}/include/${HEADER}.in"
		COMMENT "Generating header file ${HEADER}."
		VERBATIM)
endforeach(HEADER)

# Additional required packages.
find_package(LHAPDF REQUIRED)

add_library(
	sidis
	asymmetry.cpp
	bound.cpp
	constant.cpp
	cross_section.cpp
	cut.cpp
	exception.cpp
	frame.cpp
	hadronic_coeff.cpp
	kinematics.cpp
	leptonic_coeff.cpp
	math.cpp
	particle.cpp
	structure_function.cpp
	tmd.cpp
	transform.cpp
	vector.cpp
	sf_set/prokudin.cpp
	sf_set/stf.cpp
	${HEADER_LIST_SOURCE})
target_link_libraries(sidis PRIVATE MSTWPDF::mstwpdf LHAPDF Cubature::cubature)
target_include_directories(
	sidis PUBLIC
	"${Sidis_SOURCE_DIR}/include"
	"${Sidis_BINARY_DIR}/include")
target_compile_features(sidis PUBLIC cxx_std_11)
set_target_properties(
	sidis PROPERTIES
	CXX_EXTENSIONS OFF
	INTERPROCEDURAL_OPTIMIZATION ${Sidis_IPO_ENABLED})
target_compile_options(
	sidis
	PRIVATE
	$<$<OR:$<CXX_COMPILER_ID:AppleClang>,$<CXX_COMPILER_ID:Clang>,$<CXX_COMPILER_ID:GNU>>:-Wall -Wextra -pedantic>
	$<$<CXX_COMPILER_ID:MSVC>:/W4>)
add_custom_command(
	TARGET sidis POST_BUILD
	COMMAND ${CMAKE_COMMAND} -E make_directory
		"${Sidis_BINARY_DIR}/share/sidis/sf_set"
	COMMAND ${CMAKE_COMMAND} -E create_symlink
		"${Sidis_SOURCE_DIR}/share/sidis/sf_set/prokudin"
		"${Sidis_BINARY_DIR}/share/sidis/sf_set/prokudin"
	VERBATIM)

source_group(
	TREE "${Sidis_SOURCE_DIR}/include"
	PREFIX "Header files"
	FILES ${HEADER_LIST_SOURCE})

include(GNUInstallDirs)
set_target_properties(
	sidis PROPERTIES
	INSTALL_RPATH "${CMAKE_INSTALL_PREFIX}/${CMAKE_INSTALL_LIBDIR}"
	INSTALL_RPATH_USE_LINK_PATH TRUE)
target_compile_definitions(
	sidis PRIVATE
	DATADIR="${CMAKE_INSTALL_PREFIX}/${CMAKE_INSTALL_DATADIR}")
install(
	DIRECTORY "${Sidis_SOURCE_DIR}/share/sidis/"
	DESTINATION "${CMAKE_INSTALL_DATADIR}/sidis")
install(
	DIRECTORY "${Sidis_SOURCE_DIR}/include/sidis"
	DESTINATION "${CMAKE_INSTALL_INCLUDEDIR}"
	FILES_MATCHING
		PATTERN "*.hpp"
		PATTERN "*.h")
install(
	DIRECTORY "${Sidis_BINARY_DIR}/include/sidis"
	DESTINATION "${CMAKE_INSTALL_INCLUDEDIR}"
	FILES_MATCHING
		PATTERN "*.hpp"
		PATTERN "*.h")
install(
	TARGETS sidis
	LIBRARY DESTINATION "${CMAKE_INSTALL_LIBDIR}"
	ARCHIVE DESTINATION "${CMAKE_INSTALL_LIBDIR}"
	INCLUDES DESTINATION "${CMAKE_INSTALL_INCLUDEDIR}")

