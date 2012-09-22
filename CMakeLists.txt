cmake_minimum_required(VERSION 2.8)
project(MolDB)

# set the module path
set(CMAKE_MODULE_PATH ${MolDB_SOURCE_DIR}/cmake/modules)

# find and setup Boost
find_package(Boost COMPONENTS system iostreams thread REQUIRED)
include_directories(${Boost_INCLUDE_DIR})

# for boost process
include_directories(${MolDB_SOURCE_DIR}/boost)

# find and setup OpenBabel
find_package(OpenBabel2 REQUIRED)
include_directories(${OPENBABEL2_INCLUDE_DIR})

# find and setup ChemKit
find_package(ChemKit)
if(CHEMKIT_FOUND)
    include_directories(${CHEMKIT_INCLUDE_DIR})
    add_definitions(-DHAVE_CHEMKIT)
endif()

# find and setup RDKit
find_package(RDKit)
if(RDKIT_FOUND)
    include_directories(${RDKIT_INCLUDE_DIR})
    add_definitions(-DHAVE_RDKIT)
endif()

# unit tests
enable_testing()
add_subdirectory(test)

option(ENABLE_OPENMP OFF)
if(ENABLE_OPENMP)
  set(CMAKE_CXX_FLAGS "${CMAKE_CXX_FLAGS} -fopenmp")
  add_definitions(-DHAVE_OPENMP)
endif()

option(ENABLE_CUDA OFF)
if(ENABLE_CUDA)
  find_package(CUDA)
  if (CUDA_FOUND)
    include_directories(${CUDA_INCLUDE_DIR})
    add_definitions(-DHAVE_CUDA)
  endif()
endif()

set(libmoldb_srcs
    src/database.cpp
    src/messagehandler.cpp
    src/script_interface.cpp
)

add_library(moldb SHARED ${libmoldb_srcs})
if(CHEMKIT_FOUND AND RDKIT_FOUND)
    target_link_libraries(moldb ${OPENBABEL2_LIBRARIES} ${CHEMKIT_LIBRARIES} ${RDKIT_FINGERPRINTS_LIBRARY} ${RDKIT_SMILESPARSE_LIBRARY} ${Boost_LIBRARIES})
elseif(CHEMKIT_FOUND)
    target_link_libraries(moldb ${OPENBABEL2_LIBRARIES} ${CHEMKIT_LIBRARIES} ${Boost_LIBRARIES})
elseif(RDKIT_FOUND)
    target_link_libraries(moldb ${OPENBABEL2_LIBRARIES} ${RDKIT_FINGERPRINTS_LIBRARY} ${RDKIT_SMILESPARSE_LIBRARY} ${Boost_LIBRARIES})
else()
    target_link_libraries(moldb ${OPENBABEL2_LIBRARIES} ${Boost_LIBRARIES})
endif()
install(TARGETS moldb
                RUNTIME DESTINATION bin
                LIBRARY DESTINATION lib
                ARCHIVE DESTINATION lib)

add_subdirectory(tools)

#set(CPACK_SET_DESTDIR TRUE)
set(CPACK_GENERATOR "TGZ")
#set(CPACK_SOURCE_INSTALLED_DIRECTORIES ${CMAKE_SOURCE_DIR}/src;./src)
set(CPACK_SOURCE_IGNORE_FILES build;.git)
include(CPack)


