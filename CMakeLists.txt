#------------------ Caracteristicas generales del proyecto ------------------#

cmake_minimum_required(VERSION 2.8.9)

project(phoenix LANGUAGES CXX C)

set(CMAKE_MODULE_PATH "${CMAKE_SOURCE_DIR}/cmake")

set(CMAKE_BUILD_TYPE Release)

# set(CMAKE_BUILD_TYPE Debug)

set(CMAKE_BINARY_DIR ${CMAKE_SOURCE_DIR}/build)

set(EXECUTABLE_OUTPUT_PATH ${CMAKE_SOURCE_DIR}/bin)

set(LIBRARY_OUTPUT_PATH ${CMAKE_SOURCE_DIR}/lib)



# External libraries: MPI

find_package(MPI REQUIRED)


# External libraries: Boost

find_package( Boost 1.55 COMPONENTS system program_options REQUIRED )





# Bring the headers

include_directories(include ${MPI_INCLUDE_PATH})



# Compiler options for CXX

SET( CMAKE_CXX_FLAGS  "${CMAKE_CXX_FLAGS} -Wall -lm -pedantic --std=c++11" )




set ( PROJECT_LINK_LIBS
  liblatticemodel.so
  libalgebra.so
  libio.so
  libsimulation.so
  libmesh.so
  m
  # ${MPI_LIBRARIES}
  )


link_directories( ${CMAKE_SOURCE_DIR}/lib )


# # Optional: Metis decomposition

# option(WITH_METIS "Compile with METIS" ON)

# if(WITH_METIS)
#   add_definitions(-DUSE_METIS)
# endif()





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





# # Optional: CUDA

# option(WITH_CUDA "Compile with CUDA" OFF)

# if(WITH_CUDA)
  
#   find_package(CUDA QUIET REQUIRED)

#   set(CUDA_NVCC_FLAGS "${CUDA_NVCC_FLAGS} -arch=sm_30")
  
# endif()




# Directorio con fuentes

add_subdirectory (src)




# # Optional docs

# option(WITH_DOCS "Create documentation" OFF)

# if(WITH_DOCS)
#   add_subdirectory (docs)
# endif()
