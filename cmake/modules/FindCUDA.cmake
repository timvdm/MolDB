# - find CUDA
# CUDA_INCLUDE_DIR - Where to find CUDA header files (directory)
# CUDA_FINGERPRINT_LIBRARY - CUDA Fingerprint library
# CUDA_SMILESPARSE_LIBRARY - CUDA SmilesParse library
# CUDA_FOUND - Set to TRUE if we found everything (library, includes and executable)

# Copyright (c) 2010 Tim Vandermeersch
#
# Redistribution and use is allowed according to the terms of the BSD license.
# For details see the accompanying COPYING-CMAKE-SCRIPTS file.

IF(CUDA_INCLUDE_DIR AND CUDA_LIBRARIES)
    SET(CUDA_FIND_QUIETLY TRUE)
ENDIF(CUDA_INCLUDE_DIR AND CUDA_LIBRARIES)

FIND_PATH(CUDA_INCLUDE_DIR cuda_runtime.h PATHS ${CUDA_ROOT}/include /usr/local/cuda/include)

#FIND_LIBRARY(CUDA_FINGERPRINTS_LIBRARY NAMES Fingerprints PATHS ${CUDA_ROOT}/lib ${CUDA_ROOT}/build/lib)
#FIND_LIBRARY(CUDA_SMILESPARSE_LIBRARY NAMES SmilesParse PATHS ${CUDA_ROOT}/lib ${CUDA_ROOT}/build/lib)

#IF(CUDA_FINGERPRINTS_LIBRARY AND CUDA_SMILESPARSE_LIBRARY AND CUDA_INCLUDE_DIR)
IF(CUDA_INCLUDE_DIR)
    SET(CUDA_FOUND TRUE)
ENDIF()

IF(CUDA_FOUND)
    IF(NOT CUDA_FIND_QUIETLY)
        MESSAGE(STATUS "Found CUDA header file in ${CUDA_INCLUDE_DIR}")
        MESSAGE(STATUS "Found CUDA libraries: ${CUDA_FINGERPRINTS_LIBRARY} ${CUDA_SMILESPARSE_LIBRARY}")
    ENDIF(NOT CUDA_FIND_QUIETLY)
ELSE()
    IF(CUDA_FIND_REQUIRED)
      MESSAGE(FATAL_ERROR "Could not find CUDA")
    ELSE()
      MESSAGE(STATUS "Optional package CUDA was not found")
    ENDIF()
ENDIF()
