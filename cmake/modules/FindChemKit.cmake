# - find ChemKit
# CHEMKIT_INCLUDE_DIR - Where to find ChemKit header files (directory)
# CHEMKIT_LIBRARIES - ChemKit libraries
# CHEMKIT_FOUND - Set to TRUE if we found everything (library, includes and executable)

# Copyright (c) 2010 Tim Vandermeersch
#
# Redistribution and use is allowed according to the terms of the BSD license.
# For details see the accompanying COPYING-CMAKE-SCRIPTS file.

IF( CHEMKIT_INCLUDE_DIR AND CHEMKIT_LIBRARIES )
    SET(CHEMKIT_FIND_QUIETLY TRUE)
ENDIF( CHEMKIT_INCLUDE_DIR AND CHEMKIT_LIBRARIES )

FIND_PATH( CHEMKIT_INCLUDE_DIR chemkit/graphicsmoleculewireframeitem.h )

FIND_LIBRARY(CHEMKIT_LIBRARIES NAMES chemkit )

IF( CHEMKIT_LIBRARIES AND CHEMKIT_INCLUDE_DIR )
    SET( CHEMKIT_FOUND TRUE )
ENDIF( CHEMKIT_LIBRARIES AND CHEMKIT_INCLUDE_DIR )

IF( CHEMKIT_FOUND )
    IF( NOT CHEMKIT_FIND_QUIETLY )
        MESSAGE( STATUS "Found ChemKit header file in ${CHEMKIT_INCLUDE_DIR}")
        MESSAGE( STATUS "Found ChemKit libraries: ${CHEMKIT_LIBRARIES}")
    ENDIF( NOT CHEMKIT_FIND_QUIETLY )
ELSE( CHEMKIT_FOUND )
    IF( ChemKit_FIND_REQUIRED )
      MESSAGE( FATAL_ERROR "Could not find ChemKit" )
    ELSE( ChemKit_FIND_REQUIRED )
      MESSAGE( STATUS "Optional package ChemKit was not found" )
    ENDIF( ChemKit_FIND_REQUIRED )
ENDIF( CHEMKIT_FOUND )
