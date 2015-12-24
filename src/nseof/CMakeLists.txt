set(root_src Configuration.cpp DataStructures.cpp GlobalBoundaryFactory.cpp
  FlowField.cpp LinearSolver.cpp Meshsize.cpp MemoizedMesh.cpp SimpleTimer.cpp
  MultiTimer.cpp)
set(NSEOF_FLOWMODELS_SRCS flowmodels/algebraic/FlowField.cpp
  flowmodels/algebraic/FGHStencil.cpp flowmodels/algebraic/NutStencil.cpp
  flowmodels/algebraic/HStencil.cpp flowmodels/algebraic/MinimumNutStencil.cpp)
set(parallel_src parallelManagers/PetscParallelConfiguration.cpp)
set(solvers_src solvers/PetscSolver.cpp solvers/SORSolver.cpp)
set(stencils_src stencils/BFInputStencils.cpp stencils/BFStepInitStencil.cpp
  stencils/FGHStencil.cpp stencils/MaxUStencil.cpp
  stencils/MovingWallStencils.cpp stencils/NeumannBoundaryStencils.cpp
  stencils/ObstacleStencil.cpp stencils/PeriodicBoundaryStencils.cpp
  stencils/RHSStencil.cpp stencils/VelocityStencil.cpp
  walldistance/WallDistanceManager.cpp)

set(NSEOF_SRCS ${root_src} ${NSEOF_FLOWMODELS_SRCS} ${parallel_src}
  ${solvers_src} ${stencils_src})

add_library(nseof ${NSEOF_SRCS})
target_link_libraries(nseof tinyxml2 vtk ${PETSC_LIBRARIES} ANN)
