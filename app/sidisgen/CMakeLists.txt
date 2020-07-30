add_executable(sidisgen main.cpp)
target_compile_features(sidisgen PRIVATE cxx_std_11)
set_target_properties(sidisgen PROPERTIES CXX_EXTENSIONS OFF)
target_link_libraries(sidisgen PRIVATE sidis)

