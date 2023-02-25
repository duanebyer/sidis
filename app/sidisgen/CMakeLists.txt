find_package(ROOT 6.16 REQUIRED COMPONENTS Core Physics Tree CONFIG)
add_executable(sidisgen
	main.cpp
	exception.hpp
	generator.hpp generator.ipp generator.cpp
	params.hpp params.cpp
	params_format.hpp params_format.cpp
	terminal.hpp terminal.cpp
	utility.hpp)
target_link_libraries(sidisgen PRIVATE
	sidis
	Bubble::bubble
	ROOT::Core ROOT::Physics ROOT::Tree)
if(Sidis_OPENMP_ENABLED)
	if(CMAKE_VERSION VERSION_LESS 3.10)
		find_package(OpenMP 3.0 REQUIRED)
	else()
		find_package(OpenMP 3.0 REQUIRED COMPONENTS CXX)
	endif()
	if(CMAKE_VERSION VERSION_LESS 3.9)
		find_package(Threads REQUIRED)
		target_link_libraries(sidisgen PRIVATE ${OpenMP_CXX_FLAGS} Threads::Threads)
		target_compile_options(sidisgen PRIVATE ${OpenMP_CXX_FLAGS})
	else()
		target_link_libraries(sidisgen PRIVATE OpenMP::OpenMP_CXX)
	endif()
	target_compile_definitions(sidisgen PRIVATE BUBBLE_USE_OPENMP)
endif()
target_compile_features(sidisgen PRIVATE cxx_std_11)
set_target_properties(
	sidisgen PROPERTIES
	CXX_EXTENSIONS OFF
	INTERPROCEDURAL_OPTIMIZATION ${Sidis_ENABLE_IPO})
target_compile_options(
	sidisgen
	PRIVATE
	${Sidis_COMPILER_WARNINGS})

include(GNUInstallDirs)
set_target_properties(
	sidisgen PROPERTIES
	INSTALL_RPATH "${CMAKE_INSTALL_PREFIX}/${CMAKE_INSTALL_LIBDIR}"
	INSTALL_RPATH_USE_LINK_PATH TRUE)
install(
	TARGETS sidisgen
	RUNTIME DESTINATION "${CMAKE_INSTALL_BINDIR}")

