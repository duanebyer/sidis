enable_testing()

find_package(Catch2 REQUIRED)
include(CTest)
include(Catch)

add_executable(
	sidistest
	main.cpp
	test_dilog.cpp)
target_compile_features(sidistest PRIVATE cxx_std_11)
set_target_properties(sidistest PROPERTIES CXX_EXTENSIONS OFF)
target_link_libraries(sidistest PRIVATE sidis sfww Catch2::Catch2)

catch_discover_tests(sidistest)

