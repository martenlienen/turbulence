cmake_minimum_required(VERSION 2.8)

# Disable tests with -Dtests=OFF
option(tests "Build tests" ON)

option(CLUSTER "Are we building on the cluster?" OFF)

option(SCOREP "Do you want to use scorep profiling?" OFF)

option(PROFILING "Compile in extra profiling information for gprof" OFF)

if (NOT CMAKE_BUILD_TYPE)
  message(STATUS "No build type selected, default to Debug")
  set(CMAKE_BUILD_TYPE "Debug")
endif()

# Use mpic++ over g++/c++
#
# You have to specify it before the `project` call. Otherwise cmake will go into
# an infinite reconfiguration loop.
#
# The original Makefile said to use `mpiCC` on the MAC-cluster.
if (CLUSTER)
  if(SCOREP)
	get_filename_component(SCOREP "./scripts/cluster/scorep" REALPATH)
	set(CMAKE_CXX_COMPILER ${SCOREP})
   else()
    set(CMAKE_CXX_COMPILER "mpiCC")
   endif()
else()
  set(CMAKE_CXX_COMPILER "mpic++")
endif()

project("turbulence")

add_compile_options(-std=c++14)

set(CMAKE_CXX_FLAGS_DEBUG "${CMAKE_CXX_FLAGS_DEBUG} -Wall -Werror -Wpedantic")

if (PROFILING)
  set(CMAKE_CXX_FLAGS "${CMAKE_CXX_FLAGS} -pg")
  set(CMAKE_EXE_LINKER_FLAGS "${CMAKE_EXE_LINKER_FLAGS} -pg")
endif()

list(APPEND CMAKE_MODULE_PATH "./cmake/modules")

find_package(PETSc REQUIRED)
include_directories(${PETSC_INCLUDES})
add_definitions(${PETSC_DEFINITIONS})

add_subdirectory("3rdparty/tinyxml2")
include_directories("3rdparty/tinyxml2")

add_subdirectory("3rdparty/ann_1.1.2")
include_directories("3rdparty/ann_1.1.2/include")

add_subdirectory("src/vtk")
include_directories("src")

add_subdirectory("src/nseof")

add_executable(ns src/main.cpp)
target_link_libraries(ns nseof)

find_package(GTest)

if (tests AND GTEST_FOUND)
  set(TEST_SRCS test/stencils/StencilFunctionsTest.cpp)

  enable_testing()
  include_directories(${GTEST_INCLUDE_DIRS})

  add_executable(tests ${TEST_SRCS})
  target_link_libraries(tests nseof ${GTEST_BOTH_LIBRARIES})

  add_test(tests tests)
endif()
