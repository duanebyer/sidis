source_group(
	TREE "${Sidis_SOURCE_DIR}/include"
	PREFIX "Header files"
	FILES ${HEADER_LIST})

set(
	HEADER_LIST
	"${Sidis_SOURCE_DIR}/include/sidis/sidis.hpp"
	"${Sidis_SOURCE_DIR}/include/sidis/constant.hpp"
	"${Sidis_SOURCE_DIR}/include/sidis/cross_section.hpp"
	"${Sidis_SOURCE_DIR}/include/sidis/hadronic_coeff.hpp"
	"${Sidis_SOURCE_DIR}/include/sidis/kinematics.hpp"
	"${Sidis_SOURCE_DIR}/include/sidis/leptonic_coeff.hpp"
	"${Sidis_SOURCE_DIR}/include/sidis/math.hpp"
	"${Sidis_SOURCE_DIR}/include/sidis/numeric.hpp"
	"${Sidis_SOURCE_DIR}/include/sidis/vector.hpp")
add_library(
	sidis
	cross_section.cpp
	hadronic_coeff.cpp
	kinematics.cpp
	leptonic_coeff.cpp
	math.cpp
	vector.cpp
	${HEADER_LIST})
target_include_directories(sidis PUBLIC "${Sidis_SOURCE_DIR}/include")
target_compile_features(sidis PUBLIC cxx_std_11)
set_target_properties(sidis PROPERTIES CXX_EXTENSIONS OFF)
target_compile_options(
	sidis
	PRIVATE
	$<$<OR:$<CXX_COMPILER_ID:AppleClang>,$<CXX_COMPILER_ID:Clang>,$<CXX_COMPILER_ID:GNU>>:-Wall -Wextra -pedantic>
	$<$<CXX_COMPILER_ID:MSVC>:/W4>)

find_package(Mathematica REQUIRED COMPONENTS WSTP)
set(
	MODEL_WW_HEADER_LIST
	"${Sidis_SOURCE_DIR}/include/sidis/structure_function.hpp"
	"${Sidis_SOURCE_DIR}/include/sidis/sf_model/ww.hpp")
add_library(
	sfww
	sf_model/ww.cpp
	${MODEL_WW_HEADER_LIST})
target_include_directories(
	sfww
	PUBLIC
	"${Sidis_SOURCE_DIR}/include"
	PRIVATE
	"${Mathematica_WSTP_INCLUDE_DIR}")
target_link_libraries(
	sfww
	PUBLIC
	"${Mathematica_WSTP_LIBRARIES}")
if(Mathematica_WSTP_LINKER_FLAGS)
	set_target_properties(
		sfww PROPERTIES
		LINK_FLAGS "${Mathematica_WSTP_LINKER_FLAGS}")
endif()
add_custom_command(
	TARGET sfww POST_BUILD
	COMMAND ${CMAKE_COMMAND} -E make_directory
		"${Sidis_BINARY_DIR}/share/sidis/sf_model"
	COMMAND ${CMAKE_COMMAND} -E create_symlink
		"${Sidis_SOURCE_DIR}/share/sidis/sf_model/wwsidis"
		"${Sidis_BINARY_DIR}/share/sidis/sf_model/wwsidis")
target_compile_features(sfww PUBLIC cxx_std_11)
set_target_properties(sfww PROPERTIES CXX_EXTENSIONS OFF)
target_compile_options(
	sfww
	PRIVATE
	$<$<OR:$<CXX_COMPILER_ID:AppleClang>,$<CXX_COMPILER_ID:Clang>,$<CXX_COMPILER_ID:GNU>>:-Wall -Wextra -pedantic>
	$<$<CXX_COMPILER_ID:MSVC>:/W4>)

include(GNUInstallDirs)
set_target_properties(
	sidis PROPERTIES
	INSTALL_RPATH "${CMAKE_INSTALL_PREFIX}/${CMAKE_INSTALL_LIBDIR}"
	INSTALL_RPATH_USE_LINK_PATH TRUE)
set_target_properties(
	sfww PROPERTIES
	INSTALL_RPATH "${CMAKE_INSTALL_PREFIX}/${CMAKE_INSTALL_LIBDIR}"
	INSTALL_RPATH_USE_LINK_PATH TRUE)
target_compile_definitions(
	sfww PRIVATE
	DATADIR="${CMAKE_INSTALL_PREFIX}/${CMAKE_INSTALL_DATADIR}")
install(
	DIRECTORY "${Sidis_SOURCE_DIR}/share/sidis/sf_model/wwsidis"
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
	TARGETS sidis sfww
	LIBRARY DESTINATION "${CMAKE_INSTALL_LIBDIR}"
	ARCHIVE DESTINATION "${CMAKE_INSTALL_LIBDIR}"
	INCLUDES DESTINATION "${CMAKE_INSTALL_INCLUDEDIR}")

