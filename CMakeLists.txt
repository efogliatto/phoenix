#------------------ Caracteristicas generales del proyecto ------------------#

cmake_minimum_required(VERSION 2.8.9)

project(phoenix LANGUAGES CXX)

set(CMAKE_MODULE_PATH "${CMAKE_SOURCE_DIR}/cmake")

set(CMAKE_BUILD_TYPE Release)

# set(CMAKE_BUILD_TYPE Debug)

set(CMAKE_BINARY_DIR ${CMAKE_SOURCE_DIR}/build)

set(EXECUTABLE_OUTPUT_PATH ${CMAKE_SOURCE_DIR}/bin)

set(LIBRARY_OUTPUT_PATH ${CMAKE_SOURCE_DIR}/lib)


# # Bibliotecas externas: MPI

# find_package(MPI REQUIRED)



# Bring the headers

include_directories(include ${MPI_INCLUDE_PATH})



# Compiler options for CXX

SET( CMAKE_CXX_FLAGS  "${CMAKE_CXX_FLAGS} -Wall -lm -pedantic --std=c++11" )




# set ( PROJECT_LINK_LIBS
#   libio.so
#   libdictio.so
#   libbasic.so
#   liblatticemodel.so
#   libgenerallbe.so
#   libpseudopot.so
#   libpseudopot.so
#   liblimrtmodel.so
#   liblitempmrtmodel.so
#   libmymrtmodel.so
#   libxumrtmodel.so
#   libboundary.so
#   libeos.so
#   libfinitedifference.so
#   libpvtu.so
#   libasciiraw.so
#   libbinaryraw.so
#   libensight.so
#   libd2q9.so
#   libd3q15.so
#   m
#   ${MPI_LIBRARIES}
#   )




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

