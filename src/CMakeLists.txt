set(
	HEADER_LIST
	"${Sidis_SOURCE_DIR}/include/sidis/sidis.hpp")
add_library(
	sidis
	sidis.cpp
	${HEADER_LIST})
target_include_directories(sidis PUBLIC ${Sidis_SOURCE_DIR}/include/)
target_compile_features(sidis PUBLIC cxx_std_11)
set_target_properties(sidis PROPERTIES CXX_EXTENSIONS OFF)
source_group(
	TREE "${Sidis_SOURCE_DIR}/include/"
	PREFIX "Header files"
	FILES ${HEADER_LIST})

