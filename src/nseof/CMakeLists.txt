set(root_src Configuration.cpp DataStructures.cpp GlobalBoundaryFactory.cpp
  FlowField.cpp LinearSolver.cpp Meshsize.cpp MemoizedMesh.cpp SimpleTimer.cpp
  MultiTimer.cpp)
set(NSEOF_FLOWMODELS_SRCS flowmodels/turbulent/FlowField.cpp
  flowmodels/turbulent/HStencil.cpp flowmodels/turbulent/MaximumNutStencil.cpp
  flowmodels/turbulent/FGHStencil.cpp flowmodels/algebraic/NutStencil.cpp
  flowmodels/ke/BFInputStencils.cpp flowmodels/ke/FlowField.cpp
  flowmodels/ke/GlobalBoundaryFactory.cpp flowmodels/ke/InitStencil.cpp
  flowmodels/ke/KEStencil.cpp flowmodels/ke/KEStencilF.cpp
  flowmodels/ke/NutStencil.cpp flowmodels/ke/ObstacleStencil.cpp
  flowmodels/ke/RHSStencil.cpp)
set(parallel_src parallelManagers/PetscParallelConfiguration.cpp)
set(solvers_src solvers/PetscSolver.cpp solvers/SORSolver.cpp)
set(stencils_src stencils/BFInputStencils.cpp stencils/BFStepInitStencil.cpp
  stencils/FGHStencil.cpp stencils/MaxUStencil.cpp
  stencils/MovingWallStencils.cpp stencils/NeumannBoundaryStencils.cpp
  stencils/ObstacleStencil.cpp stencils/PeriodicBoundaryStencils.cpp
  stencils/RHSStencil.cpp stencils/VelocityStencil.cpp
  walldistance/WallDistanceManager.cpp geometry/GeometryManager.cpp)

set(NSEOF_SRCS ${root_src} ${NSEOF_FLOWMODELS_SRCS} ${parallel_src}
  ${solvers_src} ${stencils_src})

add_library(nseof ${NSEOF_SRCS})
target_link_libraries(nseof tinyxml2 vtk ${PETSC_LIBRARIES} ANN)
