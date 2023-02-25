set(
	HEADER_LIST_BASE
	"sidis/sidis.hpp"
	"sidis/asymmetry.hpp"
	"sidis/bound.hpp"
	"sidis/constant.hpp"
	"sidis/cross_section.hpp"
	"sidis/cut.hpp"
	"sidis/frame.hpp"
	"sidis/hadronic_coeff.hpp"
	"sidis/integ_params.hpp"
	"sidis/kinematics.hpp"
	"sidis/leptonic_coeff.hpp"
	"sidis/particle.hpp"
	"sidis/phenom.hpp"
	"sidis/structure_function.hpp"
	"sidis/tmd.hpp"
	"sidis/transform.hpp"
	"sidis/vector.hpp"
	"sidis/sf_set/mask.hpp"
	"sidis/sf_set/prokudin.hpp"
	"sidis/sf_set/test.hpp"
	"sidis/extra/exception.hpp"
	"sidis/extra/integrate.hpp"
	"sidis/extra/interpolate.hpp"
	"sidis/extra/interpolate.ipp"
	"sidis/extra/math.hpp"
	"sidis/extra/map.hpp")
set(
	HEADER_LIST_GENERATE
	"sidis/numeric.hpp"
	"sidis/version.hpp")

set(HEADER_LIST_SOURCE)
set(HEADER_LIST_DIST)
foreach(HEADER ${HEADER_LIST_BASE})
	list(APPEND HEADER_LIST_SOURCE "${Sidis_SOURCE_DIR}/include/${HEADER}")
	list(APPEND HEADER_LIST_DIST "${Sidis_SOURCE_DIR}/include/${HEADER}")
endforeach(HEADER)
foreach(HEADER ${HEADER_LIST_GENERATE})
	configure_file(
		"${Sidis_SOURCE_DIR}/include/${HEADER}.in"
		"${Sidis_BINARY_DIR}/include/${HEADER}")
	list(APPEND HEADER_LIST_SOURCE "${Sidis_SOURCE_DIR}/include/${HEADER}.in")
	list(APPEND HEADER_LIST_DIST "${Sidis_BINARY_DIR}/include/${HEADER}")
endforeach(HEADER)

find_package(GSL REQUIRED)

add_library(
	sidis
	asymmetry.cpp
	bound.cpp
	cross_section.cpp
	cut.cpp
	exception.cpp
	frame.cpp
	hadronic_coeff.cpp
	kinematics.cpp
	leptonic_coeff.cpp
	math.cpp
	particle.cpp
	phenom.cpp
	structure_function.cpp
	tmd.cpp
	transform.cpp
	vector.cpp
	sf_set/prokudin.cpp
	${HEADER_LIST_SOURCE})
target_link_libraries(sidis PRIVATE MSTWPDF::mstwpdf Cubature::cubature GSL::gsl)
target_include_directories(
	sidis PUBLIC
	"${Sidis_SOURCE_DIR}/include"
	"${Sidis_BINARY_DIR}/include")
target_compile_features(sidis PUBLIC cxx_std_11)
set_target_properties(
	sidis PROPERTIES
	CXX_EXTENSIONS OFF
	INTERPROCEDURAL_OPTIMIZATION ${Sidis_ENABLE_IPO})
target_compile_options(
	sidis
	PRIVATE
	${Sidis_COMPILER_WARNINGS})
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
get_target_property(MSTWPDF_TYPE mstwpdf TYPE)
if(MSTWPDF_TYPE STREQUAL "SHARED_LIBRARY")
	if (CMAKE_VERSION VERSION_GREATER_EQUAL 3.21)
		# TODO: As of writing, there is no way to use the MSTWPDF alias in this
		# command with CMake, even though that would be ideal.
		install(
			IMPORTED_RUNTIME_ARTIFACTS mstwpdf
			LIBRARY DESTINATION "${CMAKE_INSTALL_LIBDIR}")
	else()
		message(
			WARNING
			"MSTWPDF must be installed as a shared library manually")
	endif()
endif()

