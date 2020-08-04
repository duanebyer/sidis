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
target_include_directories(sidis PUBLIC ${Sidis_SOURCE_DIR}/include/)
target_compile_features(sidis PUBLIC cxx_std_11)
set_target_properties(sidis PROPERTIES CXX_EXTENSIONS OFF)
source_group(
	TREE "${Sidis_SOURCE_DIR}/include/"
	PREFIX "Header files"
	FILES ${HEADER_LIST})

