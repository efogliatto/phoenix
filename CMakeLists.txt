#------------------ Caracteristicas generales del proyecto ------------------#

cmake_minimum_required(VERSION 2.8.9)

project(phoenix LANGUAGES CXX C)

set(CMAKE_MODULE_PATH "${CMAKE_SOURCE_DIR}/cmake/Modules")

set(CMAKE_BUILD_TYPE Release)

set(CMAKE_BINARY_DIR ${CMAKE_SOURCE_DIR}/build)

set(EXECUTABLE_OUTPUT_PATH ${CMAKE_SOURCE_DIR}/bin)

set(LIBRARY_OUTPUT_PATH ${CMAKE_SOURCE_DIR}/lib)



# External libraries: MPI

find_package(MPI REQUIRED)



# Compiler options for CXX

SET( CMAKE_CXX_FLAGS  "${CMAKE_CXX_FLAGS} -Wall -lm -pedantic --std=c++11 -funroll-loops" )

if(CMAKE_BUILD_TYPE STREQUAL Debug)

  message("-- Compiling with debug options")

  SET( CMAKE_CXX_FLAGS  "${CMAKE_CXX_FLAGS} -g" )  
  
endif()




set ( PROJECT_LINK_LIBS
  liblatticemodel.so
  libalgebra.so
  libio.so
  libgeometry.so
  libsimulation.so
  libmesh.so
  liblatticefields.so
  liblbequations.so
  libtequations.so
  libforces.so
  libheatsources.so
  librelaxmodel.so
  m
  # ${MPI_LIBRARIES}
  )


link_directories( ${CMAKE_SOURCE_DIR}/lib )




# Optional: Metis decomposition

option(WITH_METIS "Compile with METIS" ON)

if(WITH_METIS)

  message("-- Compile with Metis")

  
  find_package( METIS REQUIRED )
  # find_package( METIS REQUIRED PATHS /users/fogliate/Programas/metis-5.1.0 )  

  include_directories(${METIS_INCLUDE_DIRS})
  
  set(LIBS ${LIBS} ${METIS_LIBRARIES})
  

  add_definitions(-DUSE_METIS)

  
else()

  message("-- Do not compile with Metis")
  
endif()



# # Optional external libraries: Boost

# option(EXT_BOOST "Compile with external BOOST" OFF)

# if(EXT_BOOST)

#   set(BOOST_ROOT "/home/fogliate/Programas/boost_1_67_0")

#   set(BOOST_INCLUDEDIR "/home/fogliate/Programas/boost_1_67_0")

# else()

#   find_package( Boost 1.55 COMPONENTS system program_options REQUIRED )
  
# endif()


# Optional external libraries: Boost

set(BOOST_BASE_DIR "system" CACHE STRING "Boost installation type")

if(BOOST_BASE_DIR STREQUAL "system")
  
  find_package( Boost 1.55 COMPONENTS system program_options REQUIRED )

else()

  set(BOOST_ROOT "${BOOST_BASE_DIR}")

  set(BOOST_INCLUDEDIR "${BOOST_BASE_DIR}")
  
  message("-- Use Boost from ${BOOST_BASE_DIR}")       
  
endif()



# Bring the headers

include_directories(include ${MPI_INCLUDE_PATH} ${BOOST_INCLUDEDIR})



# Optional: precision

set(PRECISION "double" CACHE STRING "Floating-point precision")

if(PRECISION STREQUAL "double")

  message("-- Double precision")

  add_definitions(-DDP)

  
  elseif(PRECISION STREQUAL "simple")

    add_definitions(-DSP)

    message("-- Simple precision")

  else()

    message( FATAL_ERROR "Precision not supported" )
  
endif()




# Optional: CGAL

set(CGAL_BASE_DIR "system" CACHE STRING "CGAL installation type")

if(CGAL_BASE_DIR STREQUAL "system")
  
  find_package(CGAL REQUIRED)

  add_definitions(-DUSE_CGAL)

  add_definitions(-DCGAL_USE_BASIC_VIEWER)

  message("-- Use system CGAL")

else()

  if(CGAL_BASE_DIR STREQUAL "none")

    message("-- Do not use CGAL")  

  else()

    set(CGAL_DIR "${CGAL_BASE_DIR}")    

    add_definitions(-DUSE_CGAL)

    add_definitions(-DCGAL_USE_BASIC_VIEWER)  

    message("-- Use CGAL from ${CGAL_DIR}")    
    
  endif()
  
endif()




# Sources

add_subdirectory (src)



# Optional: documentation

option(BUILD_DOC "Build documentation" ON)


# Only in Release versions

if (CMAKE_BUILD_TYPE MATCHES "^[Rr]elease")

  
  # check if Doxygen is installed

  find_package(Doxygen)

  if (DOXYGEN_FOUND)

    # set input and output files
    set(DOXYGEN_IN ${CMAKE_SOURCE_DIR}/docs/Doxyfile.in)
    set(DOXYGEN_OUT ${CMAKE_CURRENT_BINARY_DIR}/Doxyfile)

    # request to configure the file
    configure_file(${DOXYGEN_IN} ${DOXYGEN_OUT} @ONLY)
    message("Doxygen build started")
    
    # note the option ALL which allows to build the docs together with the application
    add_custom_target( docs ALL
      COMMAND ${DOXYGEN_EXECUTABLE} ${DOXYGEN_OUT}
      WORKING_DIRECTORY ${CMAKE_CURRENT_BINARY_DIR}
      COMMENT "Generating API documentation with Doxygen"
      VERBATIM )
    
  else (DOXYGEN_FOUND)
    
    message("Doxygen need to be installed to generate the doxygen documentation")
    
  endif (DOXYGEN_FOUND)
  

endif()






