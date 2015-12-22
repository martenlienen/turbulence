cmake_minimum_required(VERSION 2.8)

# Disable tests with -Dtests=OFF
option(tests "Build tests" ON)

option(CLUSTER "Are we building on the cluster?" OFF)

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
  set(CMAKE_CXX_COMPILER "mpiCC")
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

include_directories("3rdparty/ann_1.1.2/include")

# Some files expect the root directory to be on the include path
include_directories(".")

add_subdirectory("3rdparty/tinyxml2")
add_subdirectory("3rdparty/ann_1.1.2")
add_subdirectory("vtk")

set(root_src Configuration.cpp DataStructures.cpp FlowFieldLaminar.cpp
  FlowFieldTurbA.cpp GlobalBoundaryFactory.cpp LinearSolver.cpp Meshsize.cpp
  MemoizedMesh.cpp SimpleTimer.cpp MultiTimer.cpp)
set(parallel_src parallelManagers/PetscParallelConfiguration.cpp)
set(solvers_src solvers/PetscSolver.cpp solvers/SORSolver.cpp)
set(stencils_src stencils/BFInputStencils.cpp stencils/BFStepInitStencil.cpp
  stencils/FGHStencil.cpp stencils/MaxUStencil.cpp
  stencils/MovingWallStencils.cpp stencils/NeumannBoundaryStencils.cpp
  stencils/ObstacleStencil.cpp stencils/PeriodicBoundaryStencils.cpp
  stencils/FGHStencilTurb.cpp stencils/NutStencil.cpp stencils/HStencil.cpp
  stencils/MinimumNutStencil.cpp stencils/RHSStencil.cpp
  stencils/VelocityStencil.cpp walldistance/WallDistanceManager.cpp)

set(NS_EOF_SRCS ${root_src} ${parallel_src} ${solvers_src} ${stencils_src})

add_executable(ns main.cpp ${NS_EOF_SRCS})
target_link_libraries(ns tinyxml2 vtk ${PETSC_LIBRARIES} ANN)

find_package(GTest)

if (tests AND GTEST_FOUND)
  set(TEST_SRCS test/stencils/StencilFunctionsTest.cpp)

  enable_testing()
  include_directories(${GTEST_INCLUDE_DIRS})

  add_executable(tests ${TEST_SRCS} ${NS_EOF_SRCS})
  target_link_libraries(tests tinyxml2 vtk ${PETSC_LIBRARIES} ANN)

  target_link_libraries(tests ${GTEST_BOTH_LIBRARIES} ANN)

  add_test(tests tests)
endif()
