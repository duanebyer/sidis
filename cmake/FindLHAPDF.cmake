# Taken from HEP software foundation.
# - Try to find LHAPDF
# Defines:
#
#  LHAPDF_FOUND
#  LHAPDF_INCLUDE_DIR
#  LHAPDF_INCLUDE_DIRS (not cached)
#  LHAPDF_LIBRARY
#  LHAPDF_LIBRARIES (not cached)
#  LHAPDF_LIBRARY_DIRS (not cached)

find_library(LHAPDF_LIBRARY NAMES LHAPDF
             HINTS $ENV{LHAPDF_ROOT_DIR}/lib ${LHAPDF_ROOT_DIR}/lib)

find_path(LHAPDF_INCLUDE_DIR LHAPDF/LHAPDF.h
          HINTS $ENV{LHAPDF_ROOT_DIR}/include ${LHAPDF_ROOT_DIR}/include)

mark_as_advanced(LHAPDF_INCLUDE_DIR LHAPDF_LIBRARY)

# handle the QUIETLY and REQUIRED arguments and set LHAPDF_FOUND to TRUE if
# all listed variables are TRUE
include(FindPackageHandleStandardArgs)
FIND_PACKAGE_HANDLE_STANDARD_ARGS(LHAPDF DEFAULT_MSG LHAPDF_INCLUDE_DIR LHAPDF_LIBRARY)

set(LHAPDF_LIBRARIES ${LHAPDF_LIBRARY})
get_filename_component(LHAPDF_LIBRARY_DIRS ${LHAPDF_LIBRARY} PATH)

set(LHAPDF_INCLUDE_DIRS ${LHAPDF_INCLUDE_DIR})

mark_as_advanced(LHAPDF_FOUND)

