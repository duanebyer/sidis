add_executable(cross_section cross_section.cpp)
target_link_libraries(cross_section PRIVATE sidis)
target_compile_features(cross_section PRIVATE cxx_std_11)
set_target_properties(cross_section PROPERTIES CXX_EXTENSIONS OFF)
target_compile_options(
	cross_section
	PRIVATE
	$<$<OR:$<CXX_COMPILER_ID:AppleClang>,$<CXX_COMPILER_ID:Clang>,$<CXX_COMPILER_ID:GNU>>:-Wall -Wextra -pedantic>
	$<$<CXX_COMPILER_ID:MSVC>:/W4>)

