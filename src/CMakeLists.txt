set(
	HEADER_LIST_COMMON
	"${Sidis_SOURCE_DIR}/include/sidis/sidis.hpp"
	"${Sidis_SOURCE_DIR}/include/sidis/constant.hpp"
	"${Sidis_SOURCE_DIR}/include/sidis/cross_section.hpp"
	"${Sidis_SOURCE_DIR}/include/sidis/hadronic_coeff.hpp"
	"${Sidis_SOURCE_DIR}/include/sidis/integrate.hpp"
	"${Sidis_SOURCE_DIR}/include/sidis/interpolate.hpp"
	"${Sidis_SOURCE_DIR}/include/sidis/interpolate.ipp"
	"${Sidis_SOURCE_DIR}/include/sidis/kinematics.hpp"
	"${Sidis_SOURCE_DIR}/include/sidis/leptonic_coeff.hpp"
	"${Sidis_SOURCE_DIR}/include/sidis/math.hpp"
	"${Sidis_SOURCE_DIR}/include/sidis/structure_function.hpp"
	"${Sidis_SOURCE_DIR}/include/sidis/transform.hpp"
	"${Sidis_SOURCE_DIR}/include/sidis/vector.hpp"
	"${Sidis_SOURCE_DIR}/include/sidis/sf_model/ww.hpp")
set(HEADER_LIST_SOURCE ${HEADER_LIST_COMMON})
set(HEADER_LIST ${HEADER_LIST_COMMON})
list(
	APPEND HEADER_LIST_SOURCE
	"${Sidis_SOURCE_DIR}/include/sidis/numeric.hpp.in")
list(
	APPEND HEADER_LIST
	"${Sidis_SOURCE_DIR}/include/sidis/numeric.hpp")
add_library(
	sidis
	cross_section.cpp
	hadronic_coeff.cpp
	kinematics.cpp
	leptonic_coeff.cpp
	math.cpp
	transform.cpp
	vector.cpp
	sf_model/ww.cpp
	${HEADER_LIST})
target_link_libraries(sidis PRIVATE MSTWPDF::mstwpdf)
target_include_directories(
	sidis PUBLIC
	"${Sidis_SOURCE_DIR}/include"
	"${Sidis_BINARY_DIR}/include")
target_compile_features(sidis PUBLIC cxx_std_11)
set_target_properties(sidis PROPERTIES CXX_EXTENSIONS OFF)
target_compile_options(
	sidis
	PRIVATE
	$<$<OR:$<CXX_COMPILER_ID:AppleClang>,$<CXX_COMPILER_ID:Clang>,$<CXX_COMPILER_ID:GNU>>:-Wall -Wextra -pedantic>
	$<$<CXX_COMPILER_ID:MSVC>:/W4>)
add_custom_command(
	TARGET sidis POST_BUILD
	COMMAND ${CMAKE_COMMAND} -E make_directory
		"${Sidis_BINARY_DIR}/share/sidis/sf_model"
	COMMAND ${CMAKE_COMMAND} -E create_symlink
		"${Sidis_SOURCE_DIR}/share/sidis/sf_model/ww"
		"${Sidis_BINARY_DIR}/share/sidis/sf_model/ww")

source_group(
	TREE "${Sidis_SOURCE_DIR}/include"
	PREFIX "Header files"
	FILES ${HEADER_LIST_SOURCE} ${MODEL_WW_HEADER_LIST_SOURCE})

include(GNUInstallDirs)
set_target_properties(
	sidis PROPERTIES
	INSTALL_RPATH "${CMAKE_INSTALL_PREFIX}/${CMAKE_INSTALL_LIBDIR}"
	INSTALL_RPATH_USE_LINK_PATH TRUE)
target_compile_definitions(
	sidis PRIVATE
	DATADIR="${CMAKE_INSTALL_PREFIX}/${CMAKE_INSTALL_DATADIR}")
install(
	DIRECTORY "${Sidis_SOURCE_DIR}/share/sidis/sf_model/ww"
	DESTINATION "${CMAKE_INSTALL_DATADIR}/sidis/sf_model"
	FILES_MATCHING
		PATTERN ".git" EXCLUDE
		PATTERN "*.m"
		PATTERN "*.nb"
		PATTERN "*.dat")
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

