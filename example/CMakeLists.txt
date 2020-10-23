add_executable(cross_section cross_section.cpp)
target_link_libraries(cross_section PRIVATE sidis)
target_compile_features(cross_section PRIVATE cxx_std_11)
set_target_properties(cross_section PROPERTIES CXX_EXTENSIONS OFF)
target_compile_options(
	cross_section
	PRIVATE
	$<$<OR:$<CXX_COMPILER_ID:AppleClang>,$<CXX_COMPILER_ID:Clang>,$<CXX_COMPILER_ID:GNU>>:-Wall -Wextra -pedantic>
	$<$<CXX_COMPILER_ID:MSVC>:/W4>)

# ROOT is used to provide plotting for several of the examples. If not found,
# then don't build those examples.
find_package(ROOT 6.16 QUIET COMPONENTS Core Gpad Graf CONFIG)
if(NOT ROOT_FOUND)
	message(WARNING "Couldn't find ROOT, so some examples will not be built.")
endif()
if(ROOT_FOUND)
	add_executable(transverse_asymmetry transverse_asymmetry.cpp)
	target_link_libraries(transverse_asymmetry PRIVATE
		sidis
		cubature
		ROOT::Core ROOT::Gpad ROOT::Graf)
	target_compile_features(transverse_asymmetry PRIVATE cxx_std_11)
	set_target_properties(transverse_asymmetry PROPERTIES CXX_EXTENSIONS OFF)
	target_compile_options(
		transverse_asymmetry
		PRIVATE
		$<$<OR:$<CXX_COMPILER_ID:AppleClang>,$<CXX_COMPILER_ID:Clang>,$<CXX_COMPILER_ID:GNU>>:-Wall -Wextra -pedantic>
		$<$<CXX_COMPILER_ID:MSVC>:/W4>)
endif()

