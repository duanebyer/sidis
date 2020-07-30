add_executable(sidistest sidis_test.cpp)
target_compile_features(sidistest PRIVATE cxx_std_11)
set_target_properties(sidistest PROPERTIES CXX_EXTENSIONS OFF)
target_link_libraries(sidistest PRIVATE sidis)

add_test(NAME sidistest COMMAND sidistest)

