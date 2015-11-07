cmake_minimum_required(VERSION 2.8)

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

set(root_src Configuration.cpp DataStructures.cpp FlowField.cpp
  GlobalBoundaryFactory.cpp LinearSolver.cpp Meshsize.cpp SimpleTimer.cpp)
set(parallel_src parallelManagers/PetscParallelConfiguration.cpp)
set(solvers_src solvers/PetscSolver.cpp solvers/SORSolver.cpp)
set(stencils_src stencils/BFInputStencils.cpp stencils/BFStepInitStencil.cpp
  stencils/FGHStencil.cpp stencils/MaxUStencil.cpp stencils/VTKStencil.cpp
  stencils/MovingWallStencils.cpp stencils/NeumannBoundaryStencils.cpp
  stencils/ObstacleStencil.cpp stencils/PeriodicBoundaryStencils.cpp
  stencils/RHSStencil.cpp stencils/VelocityStencil.cpp)

add_executable(ns main.cpp ${root_src} ${parallel_src} ${solvers_src}
  ${stencils_src})
target_link_libraries(ns tinyxml2 ${PETSC_LIBRARIES})
