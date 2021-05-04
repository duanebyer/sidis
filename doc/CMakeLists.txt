if(Sidis_BUILD_DOXYGEN)
	find_package(Doxygen REQUIRED)
	# TODO: Use the CMake version 3.9 method if possible.
	configure_file(Doxyfile "${Sidis_BINARY_DIR}/doc/Doxyfile")
	add_custom_target(
		doc ALL
		COMMAND ${DOXYGEN_EXECUTABLE} "${Sidis_BINARY_DIR}/doc/Doxyfile"
		WORKING_DIRECTORY "${Sidis_BINARY_DIR}/doc"
		COMMENT "Generating documentation with Doxygen."
		VERBATIM)
	include(GNUInstallDirs)
	install(
		DIRECTORY "${Sidis_BINARY_DIR}/doc"
		DESTINATION "${CMAKE_INSTALL_DATADIR}/sidis/doc"
		FILES_MATCHING
			PATTERN "html/*"
			PATTERN "latex/*")
endif()

if(Sidis_BUILD_DOCS)
	find_package(LATEX REQUIRED)
	# Right now, there is no other documentation to be built.
endif()

