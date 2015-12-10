cmake_minimum_required(VERSION 2.8)

# Disable tests with -Dtests=OFF
option(tests "Build tests" ON)

# Use mpic++ over g++/c++
#
# You have to specify it before the `project` call. Otherwise cmake will go into
# an infinite reconfiguration loop.
#
# The original Makefile said to use `mpiCC` on the MAC-cluster.
set(CMAKE_CXX_COMPILER "mpic++")

project("turbulence")

add_compile_options(-Wall -Werror -Wpedantic -O3 -std=c++14)

list(APPEND CMAKE_MODULE_PATH "./cmake/modules")

find_package(PETSc REQUIRED)
include_directories(${PETSC_INCLUDES})
add_definitions(${PETSC_DEFINITIONS})

# Some files expect the root directory to be on the include path
include_directories(".")

add_subdirectory("3rdparty/tinyxml2")
add_subdirectory("vtk")

set(root_src Configuration.cpp DataStructures.cpp FlowField.cpp
  GlobalBoundaryFactory.cpp LinearSolver.cpp Meshsize.cpp SimpleTimer.cpp)
set(parallel_src parallelManagers/PetscParallelConfiguration.cpp)
set(solvers_src solvers/PetscSolver.cpp solvers/SORSolver.cpp)
set(stencils_src stencils/BFInputStencils.cpp stencils/BFStepInitStencil.cpp
  stencils/FGHStencil.cpp stencils/MaxUStencil.cpp
  stencils/MovingWallStencils.cpp stencils/NeumannBoundaryStencils.cpp
  stencils/ObstacleStencil.cpp stencils/PeriodicBoundaryStencils.cpp
  stencils/RHSStencil.cpp stencils/VelocityStencil.cpp)

set(NS_EOF_SRCS ${root_src} ${parallel_src} ${solvers_src} ${stencils_src})

add_executable(ns main.cpp ${NS_EOF_SRCS})
target_link_libraries(ns tinyxml2 vtk ${PETSC_LIBRARIES})

find_package(GTest)

if (tests AND GTEST_FOUND)
  set(TEST_SRCS)

  enable_testing()
  include_directories(${GTEST_INCLUDE_DIRS})

  add_executable(tests ${TEST_SRCS} ${NS_EOF_SRCS})
  target_link_libraries(tests tinyxml2 vtk ${PETSC_LIBRARIES})

  target_link_libraries(tests ${GTEST_BOTH_LIBRARIES})
  add_test(tests tests)
endif()
