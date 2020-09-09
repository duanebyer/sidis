add_executable(sidisgen main.cpp)
target_link_libraries(sidisgen PRIVATE sidis)
target_compile_features(sidisgen PRIVATE cxx_std_11)
set_target_properties(sidisgen PROPERTIES CXX_EXTENSIONS OFF)
target_compile_options(
	sidisgen
	PRIVATE
	$<$<OR:$<CXX_COMPILER_ID:AppleClang>,$<CXX_COMPILER_ID:Clang>,$<CXX_COMPILER_ID:GNU>>:-Wall -Wextra -pedantic>
	$<$<CXX_COMPILER_ID:MSVC>:/W4>)

include(GNUInstallDirs)
set_target_properties(
	sidisgen PROPERTIES
	INSTALL_RPATH "${CMAKE_INSTALL_PREFIX}/${CMAKE_INSTALL_LIBDIR}"
	INSTALL_RPATH_USE_LINK_PATH TRUE)
install(
	TARGETS sidisgen
	RUNTIME DESTINATION "${CMAKE_INSTALL_BINDIR}")

