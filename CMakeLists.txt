#set the minimum cmake version
cmake_minimum_required(VERSION 3.10)


#define project name
project(FVS2D)


#define language
enable_language(Fortran)


#set executable file location (build/run/)
set(CMAKE_RUNTIME_OUTPUT_DIRECTORY ${CMAKE_BINARY_DIR}/run)


#check for MPI 
find_package(MPI REQUIRED)
if (MPI_FOUND)
  include_directories(SYSTEM ${MPI_Fortran_INCLUDE_PATH})
else (MPI_FOUND)
  message(SEND_ERROR "This application cannot compile without MPI")
endif (MPI_FOUND)


#define flags before add src directory
if (${CMAKE_Fortran_COMPILER_ID} MATCHES "Intel")
  set(CMAKE_Fortran_FLAGS "-r8 -Ofast -fpp -convert big_endian -traceback -qopenmp -mcmodel=medium")
elseif (${CMAKE_Fortran_COMPILER_ID} MATCHES "Cray")
  set(CMAKE_Fortran_FLAGS "-s real64 -Ofast -e P -h byteswapio -fopenmp")
else()
  set(CMAKE_Fortran_FLAGS "-r8 -Ofast -fpp -convert big_endian -traceback -qopenmp -mcmodel=medium")
endif()


#add directory for source files
add_subdirectory(src)


#add directory for libraries
add_subdirectory(libs/BLAS)


#print some info for user
if(CMAKE_Fortran_COMPILER_LOADED)
  message(STATUS "Compiler path: ${CMAKE_Fortran_COMPILER}")
  message(STATUS "Compiler version: ${CMAKE_Fortran_COMPILER_VERSION}")
  message(STATUS "Compiler ID: ${CMAKE_Fortran_COMPILER_ID}")
  message(STATUS "Compiler flags: ${CMAKE_Fortran_FLAGS}")
endif()
