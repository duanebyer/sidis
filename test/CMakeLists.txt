enable_testing()

find_package(Catch2 2.7.0 REQUIRED)
include(Catch)

add_executable(
	sidistest
	main.cpp
	test_cross_section.cpp
	test_interpolate.cpp
	test_kinematics.cpp
	test_math.cpp
	test_vector.cpp
	phase_space_generator.cpp)
target_include_directories(sidistest PRIVATE ${Sidis_SOURCE_DIR}/test)
target_link_libraries(sidistest PRIVATE sidis Catch2::Catch2)
add_custom_command(
	TARGET sidistest POST_BUILD
	COMMAND ${CMAKE_COMMAND} -E create_symlink
		"${Sidis_SOURCE_DIR}/test/data"
		"${Sidis_BINARY_DIR}/test/data")
target_compile_features(sidistest PRIVATE cxx_std_11)
set_target_properties(
	sidistest PROPERTIES
	CXX_EXTENSIONS OFF
	INTERPROCEDURAL_OPTIMIZATION ${Sidis_ENABLE_IPO})

catch_discover_tests(sidistest)

