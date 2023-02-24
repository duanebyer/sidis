add_executable(cross_section cross_section.cpp)
target_link_libraries(cross_section PRIVATE sidis)
target_compile_features(cross_section PRIVATE cxx_std_11)
set_target_properties(
	cross_section PROPERTIES
	CXX_EXTENSIONS OFF
	INTERPROCEDURAL_OPTIMIZATION ${Sidis_ENABLE_IPO})
target_compile_options(
	cross_section
	PRIVATE
	${Sidis_COMPILER_WARNINGS})

add_executable(alpha_qed alpha_qed.cpp)
target_link_libraries(alpha_qed PRIVATE sidis)
target_compile_features(alpha_qed PRIVATE cxx_std_11)
set_target_properties(
	alpha_qed PROPERTIES
	CXX_EXTENSIONS OFF
	INTERPROCEDURAL_OPTIMIZATION ${Sidis_ENABLE_IPO})
target_compile_options(
	alpha_qed
	PRIVATE
	${Sidis_COMPILER_WARNINGS})

add_executable(random_phase_space random_phase_space.cpp)
target_link_libraries(random_phase_space PRIVATE sidis)
target_compile_features(random_phase_space PRIVATE cxx_std_11)
set_target_properties(
	random_phase_space PROPERTIES
	CXX_EXTENSIONS OFF
	INTERPROCEDURAL_OPTIMIZATION ${Sidis_ENABLE_IPO})
target_compile_options(
	random_phase_space
	PRIVATE
	${Sidis_COMPILER_WARNINGS})

add_executable(quick_start quick_start.cpp)
target_link_libraries(quick_start PRIVATE sidis)
target_compile_features(quick_start PRIVATE cxx_std_11)
set_target_properties(
	quick_start PROPERTIES
	CXX_EXTENSIONS OFF
	INTERPROCEDURAL_OPTIMIZATION ${Sidis_ENABLE_IPO})
target_compile_options(
	quick_start
	PRIVATE
	${Sidis_COMPILER_WARNINGS})

# ROOT is used to provide plotting for several of the examples. If not found,
# then don't build those examples.
find_package(ROOT 6.16 COMPONENTS Core Gpad Graf CONFIG)
if(NOT ROOT_FOUND)
	message(WARNING "Couldn't find ROOT, so some examples will not be built.")
endif()
if(ROOT_FOUND)
	add_executable(plot_cutoff plot_cutoff.cpp)
	target_link_libraries(plot_cutoff PRIVATE
		sidis
		ROOT::Core ROOT::Gpad ROOT::Graf)
	target_compile_features(plot_cutoff PRIVATE cxx_std_11)
	set_target_properties(
		plot_cutoff PROPERTIES
		CXX_EXTENSIONS OFF
		INTERPROCEDURAL_OPTIMIZATION ${Sidis_ENABLE_IPO})
	target_compile_options(
		plot_cutoff
		PRIVATE
		${Sidis_COMPILER_WARNINGS})

	add_executable(plot_ssa plot_ssa.cpp)
	target_link_libraries(plot_ssa PRIVATE
		sidis
		ROOT::Core ROOT::Gpad ROOT::Graf)
	target_compile_features(plot_ssa PRIVATE cxx_std_11)
	set_target_properties(
		plot_ssa PROPERTIES
		CXX_EXTENSIONS OFF
		INTERPROCEDURAL_OPTIMIZATION ${Sidis_ENABLE_IPO})
	target_compile_options(
		plot_ssa
		PRIVATE
		${Sidis_COMPILER_WARNINGS})

	add_executable(plot_phase_space plot_phase_space.cpp)
	target_link_libraries(plot_phase_space PRIVATE
		sidis
		ROOT::Core ROOT::Gpad ROOT::Graf)
	target_compile_features(plot_phase_space PRIVATE cxx_std_11)
	set_target_properties(
		plot_phase_space PROPERTIES
		CXX_EXTENSIONS OFF
		INTERPROCEDURAL_OPTIMIZATION ${Sidis_ENABLE_IPO})
	target_compile_options(
		plot_phase_space
		PRIVATE
		${Sidis_COMPILER_WARNINGS})

	add_executable(process_events process_events.cpp)
	target_link_libraries(process_events PRIVATE
		ROOT::Core ROOT::Gpad ROOT::Graf ROOT::Tree)
	target_compile_features(process_events PRIVATE cxx_std_11)
	set_target_properties(
		process_events PROPERTIES
		CXX_EXTENSIONS OFF
		INTERPROCEDURAL_OPTIMIZATION ${Sidis_ENABLE_IPO})
	target_compile_options(
		process_events
		PRIVATE
		${Sidis_COMPILER_WARNINGS})
endif()

