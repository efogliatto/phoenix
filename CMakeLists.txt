#------------------ Caracteristicas generales del proyecto ------------------#

cmake_minimum_required(VERSION 2.8.9)

project(phoenix LANGUAGES CXX C)

set(CMAKE_MODULE_PATH "${CMAKE_SOURCE_DIR}/cmake")

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



# Modules path

set(CMAKE_MODULE_PATH ${CMAKE_MODULE_PATH} "${CMAKE_SOURCE_DIR}/cmake/Modules/")


# Optional: Metis decomposition

option(WITH_METIS "Compile with METIS" ON)

if(WITH_METIS)

  message("-- Compile with Metis")

  find_package( METIS REQUIRED )

  include_directories(${METIS_INCLUDE_DIRS})
  
  set(LIBS ${LIBS} ${METIS_LIBRARIES})

  add_definitions(-DUSE_METIS)
  
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







# Directorio con fuentes

add_subdirectory (src)


