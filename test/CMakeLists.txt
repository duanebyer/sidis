enable_testing()

find_package(Catch2 REQUIRED)
include(CTest)
include(Catch)

add_executable(
	sidistest
	main.cpp
	test_math.cpp
	test_cross_section.cpp)
target_include_directories(sidistest PRIVATE ${Sidis_SOURCE_DIR}/test)
target_link_libraries(sidistest PRIVATE sidis sfww Catch2::Catch2)
add_custom_command(
	TARGET sidistest POST_BUILD
	COMMAND ${CMAKE_COMMAND} -E create_symlink
		"${Sidis_SOURCE_DIR}/test/data"
		"${Sidis_BINARY_DIR}/test/data")
target_compile_features(sidistest PRIVATE cxx_std_11)
set_target_properties(sidistest PROPERTIES CXX_EXTENSIONS OFF)

catch_discover_tests(sidistest)

