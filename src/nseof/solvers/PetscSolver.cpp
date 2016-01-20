#include "PetscSolver.h"
#include "nseof/MultiTimer.h"

namespace nseof {

// This function returns the ranges to work on the pressure with the
// non-boundary stencil.
// Since the domain PETSc deals with has an additional layer of cells, the size
// is clipped to
// ignore them.
void createLimits(Parameters &parameters, DM &da, int *limitsX, int *limitsY,
                  int *limitsZ) {
  // Location of the first element and sizes of the subdomain in each dimension.
  PetscInt firstX, lengthX, firstY, lengthY, firstZ, lengthZ;

  DMDAGetCorners(da, &firstX, &firstY, &firstZ, &lengthX, &lengthY, &lengthZ);

  int rank;
  MPI_Comm_rank(PETSC_COMM_WORLD, &rank);

  // Preliminary set for the iteration domain
  limitsX[0] = firstX;
  limitsX[1] = firstX + lengthX;
  limitsY[0] = firstY;
  limitsY[1] = firstY + lengthY;
  limitsZ[0] = firstZ;
  limitsZ[1] = firstZ + lengthZ;

  // The indices are used to determine the local range of iterations. If the
  // domain is not in a
  // global boundary, we don't have to consider special conditions for the
  // boundaries. Otherwise,
  // these are shifted so that the body iterations don't touch the boundaries,
  // which are treated
  // with a differently
  // Check left wall
  if (parameters.parallel.indices[0] == 0) {
    limitsX[0]++;
  }

  // Check right wall
  if (parameters.parallel.indices[0] ==
      parameters.parallel.numProcessors[0] - 1) {
    limitsX[1]--;
  }

  // Check bottom wall
  if (parameters.parallel.indices[1] == 0) {
    limitsY[0]++;
  }

  // Check top wall
  if (parameters.parallel.indices[1] ==
      parameters.parallel.numProcessors[1] - 1) {
    limitsY[1]--;
  }

  // Check front wall
  if (parameters.parallel.indices[2] == 0) {
    limitsZ[0]++;
  }

  // Check back wall
  if (parameters.parallel.indices[2] ==
      parameters.parallel.numProcessors[2] - 1) {
    limitsZ[1]--;
  }
}

PetscUserCtx::PetscUserCtx(Parameters &parameters, FlowField &flowField)
    : _parameters(parameters), _flowField(flowField) {}

Parameters &PetscUserCtx::getParameters() { return _parameters; }

FlowField &PetscUserCtx::getFlowField() { return _flowField; }

#if (PETSC_VERSION_MAJOR == 3 && PETSC_VERSION_MINOR >= 5)
PetscErrorCode computeMatrix2D(KSP ksp, Mat A, Mat pc, void *ctx);

PetscErrorCode computeMatrix3D(KSP ksp, Mat A, Mat pc, void *ctx);
#else
PetscErrorCode computeMatrix2D(KSP ksp, Mat A, Mat pc,
                               MatStructure *matStructure, void *ctx);

PetscErrorCode computeMatrix3D(KSP ksp, Mat A, Mat pc,
                               MatStructure *matStructure, void *ctx);

#endif

PetscErrorCode computeRHS2D(KSP ksp, Vec b, void *ctx);
PetscErrorCode computeRHS3D(KSP ksp, Vec b, void *ctx);

PetscSolver::PetscSolver(FlowField &flowField, Parameters &parameters)
    : LinearSolver(flowField, parameters), _ctx(parameters, flowField) {
// Set the type of boundary nodes of the system
#if ((PETSC_VERSION_MAJOR == 3) && (PETSC_VERSION_MINOR >= 5))
  DMBoundaryType bx = DM_BOUNDARY_NONE, by = DM_BOUNDARY_NONE,
                 bz = DM_BOUNDARY_NONE;

  if (parameters.walls.typeLeft == PERIODIC) {
    bx = DM_BOUNDARY_PERIODIC;
  }
  if (parameters.walls.typeBottom == PERIODIC) {
    by = DM_BOUNDARY_PERIODIC;
  }
  if (parameters.walls.typeFront == PERIODIC) {
    bz = DM_BOUNDARY_PERIODIC;
  }

#else
  DMDABoundaryType bx = DMDA_BOUNDARY_NONE, by = DMDA_BOUNDARY_NONE,
                   bz = DMDA_BOUNDARY_NONE;

  if (parameters.walls.typeLeft == PERIODIC) {
    bx = DMDA_BOUNDARY_PERIODIC;
  }
  if (parameters.walls.typeBottom == PERIODIC) {
    by = DMDA_BOUNDARY_PERIODIC;
  }
  if (parameters.walls.typeFront == PERIODIC) {
    bz = DMDA_BOUNDARY_PERIODIC;
  }

#endif

  KSPCreate(PETSC_COMM_WORLD, &_ksp);
  PCCreate(PETSC_COMM_WORLD, &_pc);
#if (PETSC_VERSION_MAJOR == 3 && PETSC_VERSION_MINOR >= 5)
  PetscErrorCode (*computeMatrix)(KSP, Mat, Mat, void *) = NULL;
#else
  PetscErrorCode (*computeMatrix)(KSP, Mat, Mat, MatStructure *, void *) = NULL;
#endif
  if (_parameters.geometry.dim == 2) {
    computeMatrix = computeMatrix2D;
    DMDACreate2d(PETSC_COMM_WORLD, bx, by, DMDA_STENCIL_STAR,
                 parameters.geometry.sizeX + 2, parameters.geometry.sizeY + 2,
                 parameters.parallel.numProcessors[0],
                 parameters.parallel.numProcessors[1], 1, 2,
                 _parameters.parallel.sizes[0], _parameters.parallel.sizes[1],
                 &_da);
  } else if (_parameters.geometry.dim == 3) {
    computeMatrix = computeMatrix3D;
    DMDACreate3d(PETSC_COMM_WORLD, bx, by, bz, DMDA_STENCIL_STAR,
                 parameters.geometry.sizeX + 2, parameters.geometry.sizeY + 2,
                 parameters.geometry.sizeZ + 2,
                 parameters.parallel.numProcessors[0],
                 parameters.parallel.numProcessors[1],
                 parameters.parallel.numProcessors[2], 1,
                 2,  // Degrees of freedom and stencil length
                 _parameters.parallel.sizes[0], _parameters.parallel.sizes[1],
                 _parameters.parallel.sizes[2], &_da);
  }

  // Find out what are the corners of the subdomain
  DMDAGetCorners(_da, &_firstX, &_firstY, &_firstZ, &_lengthX, &_lengthY,
                 &_lengthZ);

  // Current function to assign the limits
  createLimits(parameters, _da, _limitsX, _limitsY, _limitsZ);

  // Set the rank in the context. Necessary since PETSc declares the rank of the
  // neighbors under
  // periodic conditions as the rank of the process itself. So the rank must be
  // known to properly
  // set matrices and RHS vectors
  int rank;
  MPI_Comm_rank(PETSC_COMM_WORLD, &rank);

  // Set offsets to fix where the results of the pressure will be written in the
  // flow field
  if (_firstX == 0) {
    _offsetX = 1;
  } else {
    _offsetX = 2;
  }

  if (_firstY == 0) {
    _offsetY = 1;
  } else {
    _offsetY = 2;
  }

  if (_firstZ == 0) {
    _offsetZ = 1;
  } else {
    _offsetZ = 2;
  }

  // Set a pointer to the limits in the context, so that they can be used by the
  // function
  _ctx.setLimits(_limitsX, _limitsY, _limitsZ);

  // Determine whether a process writes a boundary on the system.
  // Right now, it only depends on the position of the array. The identity of
  // the neighbors will
  // become significant only in the iterators, where it will apply a boundary
  // condition only if
  // the rank is invalid, so that there is no neighbor
  _ctx.setAsBoundary = 0;
  if (parameters.parallel.indices[0] == 0) {
    _ctx.setAsBoundary += LEFT_WALL_BIT;
  }
  if (parameters.parallel.indices[0] ==
      parameters.parallel.numProcessors[0] - 1) {
    _ctx.setAsBoundary += RIGHT_WALL_BIT;
  }
  if (parameters.parallel.indices[1] == 0) {
    _ctx.setAsBoundary += BOTTOM_WALL_BIT;
  }
  if (parameters.parallel.indices[1] ==
      parameters.parallel.numProcessors[1] - 1) {
    _ctx.setAsBoundary += TOP_WALL_BIT;
  }
  if (parameters.parallel.indices[2] == 0) {
    _ctx.setAsBoundary += FRONT_WALL_BIT;
  }
  if (parameters.parallel.indices[2] ==
      parameters.parallel.numProcessors[2] - 1) {
    _ctx.setAsBoundary += BACK_WALL_BIT;
  }

  // Set displacements to deal with periodic boundaries if necessary.
  // If the boundary is periodic, it will take information from positions beyond
  // the ghost cells,
  // since they are used only for parallel communication. Otherwise, PETSc deals
  // with
  // communication of the pressure.
  int Nx = parameters.geometry.sizeX + 2;
  int Ny = parameters.geometry.sizeY + 2;
  int Nz = parameters.geometry.sizeZ + 2;

  if (parameters.walls.typeLeft == PERIODIC) {
    _ctx.displacement[0] = -2;
    _ctx.displacement[1] = Nx + 1;
  } else {
    _ctx.displacement[0] = 1;
    _ctx.displacement[1] = Nx - 2;
  }

  if (parameters.walls.typeBottom == PERIODIC) {
    _ctx.displacement[2] = -2;
    _ctx.displacement[3] = Ny + 1;
  } else {
    _ctx.displacement[2] = 1;
    _ctx.displacement[3] = Ny - 2;
  }

  if (parameters.walls.typeFront == PERIODIC) {
    _ctx.displacement[4] = -2;
    _ctx.displacement[5] = Nz + 1;
  } else {
    _ctx.displacement[4] = 1;
    _ctx.displacement[5] = Nz - 2;
  }

  DMCreateGlobalVector(_da, &_x);
  KSPSetDM(_ksp, _da);
  KSPSetComputeOperators(_ksp, computeMatrix, &_ctx);

  KSPSetFromOptions(_ksp);
  KSPSetInitialGuessNonzero(_ksp, PETSC_TRUE);

  PetscBool hasKsp;
  PetscBool hasPc;
  PetscOptionsHasName(NULL, "-ksp_type", &hasKsp);
  PetscOptionsHasName(NULL, "-pc_type", &hasPc);

  if (!hasKsp) {
    KSPSetType(_ksp, KSPFGMRES);
  }

  int comm_size;
  MPI_Comm_size(PETSC_COMM_WORLD, &comm_size);

  if (!hasPc) {
    if (comm_size == 1) {
      // if serial
      PCSetType(_pc, PCILU);
      PCFactorSetLevels(_pc, 1);
    } else {
      // if parallel
      PCSetType(_pc, PCASM);
    }

    KSPSetPC(_ksp, _pc);
  }

  KSPSetUp(_ksp);

  if (!hasPc && comm_size > 1) {
    KSP *subksp;
    PC subpc;

    PCASMGetSubKSP(_pc, NULL, NULL, &subksp);
    KSPGetPC(subksp[0], &subpc);

    PetscBool has_fl;
    PetscBool has_sub_type;
    PetscOptionsHasName(NULL, "-sub_pc_factor_levels", &has_fl);
    PetscOptionsHasName(NULL, "-sub_pc_type", &has_sub_type);

    if (!has_sub_type) {
      PCSetType(subpc, PCILU);
    }

    if (!has_fl) {
      PCFactorSetLevels(subpc, 1);
    }
  }
}

void PetscSolver::init() {
  // initialize PETSc with the same pressure field as in flow field

  MPI_Barrier(PETSC_COMM_WORLD);
  ScalarField &pressure = _flowField.getPressure();

  if (_parameters.geometry.dim == 2) {
    PetscScalar **array;

    // get array of vector
    DMDAVecGetArray(_da, _x, &array);

    // fill array with new values (from old simulation)
    for (int j = _firstY; j < _firstY + _lengthY; j++) {
      for (int i = _firstX; i < _firstX + _lengthX; i++) {
        array[j][i] =
            pressure.getScalar(i - _firstX + _offsetX, j - _firstY + _offsetY);
      }
    }

    // write array back to vector
    DMDAVecRestoreArray(_da, _x, &array);

  } else if (_parameters.geometry.dim == 3) {
    // same form 3D

    PetscScalar ***array;
    DMDAVecGetArray(_da, _x, &array);

    for (int k = _firstZ; k < _firstZ + _lengthZ; k++) {
      for (int j = _firstY; j < _firstY + _lengthY; j++) {
        for (int i = _firstX; i < _firstX + _lengthX; i++) {
          array[k][j][i] =
              pressure.getScalar(i - _firstX + _offsetX, j - _firstY + _offsetY,
                                 k - _firstZ + _offsetZ);
        }
      }
    }
    DMDAVecRestoreArray(_da, _x, &array);
  }
}

void PetscSolver::solve() {
  MultiTimer *timer = MultiTimer::get();
  ScalarField &pressure = _flowField.getPressure();

  if (_parameters.geometry.dim == 2) {
    timer->start("poisson");

    KSPSetComputeRHS(_ksp, computeRHS2D, &_ctx);
    KSPSetComputeOperators(_ksp, computeMatrix2D, &_ctx);
    KSPSolve(_ksp, PETSC_NULL, _x);

    timer->stop("poisson");

    // Then extract the information
    PetscScalar **array;
    DMDAVecGetArray(_da, _x, &array);

    for (int j = _firstY; j < _firstY + _lengthY; j++) {
      for (int i = _firstX; i < _firstX + _lengthX; i++) {
        pressure.getScalar(i - _firstX + _offsetX, j - _firstY + _offsetY) =
            array[j][i];
      }
    }
    DMDAVecRestoreArray(_da, _x, &array);
  } else if (_parameters.geometry.dim == 3) {
    timer->start("poisson");

    KSPSetComputeRHS(_ksp, computeRHS3D, &_ctx);
    KSPSetComputeOperators(_ksp, computeMatrix3D, &_ctx);
    KSPSolve(_ksp, PETSC_NULL, _x);

    timer->stop("poisson");

    // Then extract the information
    PetscScalar ***array;
    DMDAVecGetArray(_da, _x, &array);

    for (int k = _firstZ; k < _firstZ + _lengthZ; k++) {
      for (int j = _firstY; j < _firstY + _lengthY; j++) {
        for (int i = _firstX; i < _firstX + _lengthX; i++) {
          pressure.getScalar(i - _firstX + _offsetX, j - _firstY + _offsetY,
                             k - _firstZ + _offsetZ) = array[k][j][i];
        }
      }
    }
    DMDAVecRestoreArray(_da, _x, &array);
  }
}

#if (PETSC_VERSION_MAJOR == 3 && PETSC_VERSION_MINOR >= 5)
PetscErrorCode computeMatrix2D(KSP ksp, Mat A, Mat pc, void *ctx) {
#else
PetscErrorCode computeMatrix2D(KSP ksp, Mat A, Mat pc,
                               MatStructure *matStructure, void *ctx) {
#endif

  PetscUserCtx *context = (PetscUserCtx *)ctx;
  Parameters &parameters = context->getParameters();

  IntScalarField &flags = context->getFlowField().getFlags();

  int *limitsX, *limitsY, *limitsZ;
  context->getLimits(&limitsX, &limitsY, &limitsZ);

  PetscScalar stencilValues[5];
  MatStencil row, column[5];

  PetscInt i, j, Nx, Ny;

  Nx = parameters.geometry.sizeX + 2;
  Ny = parameters.geometry.sizeY + 2;

  std::cout << "Limits= " << limitsX[0] << ", " << limitsX[1] << "; "
            << limitsY[0] << " , " << limitsY[1] << std::endl;
  // Loop for inner nodes
  for (j = limitsY[0]; j < limitsY[1]; j++) {
    for (i = limitsX[0]; i < limitsX[1]; i++) {
      row.i = i;
      row.j = j;
      // convert matrix to Cartesian grid coordinates
      const int cellIndexX = i - limitsX[0] + 2;
      const int cellIndexY = j - limitsY[0] + 2;

      const int obstacle = flags.getValue(cellIndexX, cellIndexY);

      if ((obstacle & OBSTACLE_SELF) == 0) {  // If we have a fluid cell
        const FLOAT dx_0 = parameters.meshsize->getDx(cellIndexX, cellIndexY);
        const FLOAT dx_M1 =
            parameters.meshsize->getDx(cellIndexX - 1, cellIndexY);
        const FLOAT dx_P1 =
            parameters.meshsize->getDx(cellIndexX + 1, cellIndexY);
        const FLOAT dy_0 = parameters.meshsize->getDy(cellIndexX, cellIndexY);
        const FLOAT dy_M1 =
            parameters.meshsize->getDy(cellIndexX, cellIndexY - 1);
        const FLOAT dy_P1 =
            parameters.meshsize->getDy(cellIndexX, cellIndexY + 1);

        const FLOAT dx_L = 0.5 * (dx_0 + dx_M1);
        const FLOAT dx_R = 0.5 * (dx_0 + dx_P1);
        const FLOAT dx_Bo = 0.5 * (dy_0 + dy_M1);
        const FLOAT dx_T = 0.5 * (dy_0 + dy_P1);

        // Definition of values: set general formulation for laplace operator
        // here, based on arbitrary meshsizes
        stencilValues[1] = 2.0 / (dx_L * (dx_L + dx_R));    // left
        stencilValues[0] = 2.0 / (dx_R * (dx_L + dx_R));    // right
        stencilValues[3] = 2.0 / (dx_T * (dx_T + dx_Bo));   // top
        stencilValues[4] = 2.0 / (dx_Bo * (dx_T + dx_Bo));  // bottom
        stencilValues[2] =
            -2.0 / (dx_R * dx_L) - 2.0 / (dx_T * dx_Bo);  // center

        // Definition of positions. Order must correspond to values
        column[0].i = i + 1;
        column[0].j = j;
        column[1].i = i - 1;
        column[1].j = j;
        column[2].i = i;
        column[2].j = j;
        column[3].i = i;
        column[3].j = j + 1;
        column[4].i = i;
        column[4].j = j - 1;

        MatSetValuesStencil(A, 1, &row, 5, column, stencilValues,
                            INSERT_VALUES);
      } else if (obstacle !=
                 OBSTACLE_SELF + OBSTACLE_LEFT + OBSTACLE_RIGHT + OBSTACLE_TOP +
                     OBSTACLE_BOTTOM) {  // Not fluid, but fluid somewhere
                                         // around
        int counter = 0;  // This will contain how many neighbours are fluid
        int counter_fluid = 0;
        // std::cout<<"found obstacle cell"<<std::endl;
        // TODO here, variable meshwidth might have to be considered
        if ((obstacle & OBSTACLE_LEFT) == 0) {  // If there is fluid to the left
          stencilValues[counter] = 1.0;
          column[counter].i = i - 1;
          column[counter].j = j;
          counter++;  // We have just identified a fuid cell and prepared to
                      // average
          counter_fluid++;
        } else {
          stencilValues[counter] = 0.0;
          column[counter].i = i - 1;
          column[counter].j = j;
          counter++;
        }
        if ((obstacle & OBSTACLE_RIGHT) == 0) {
          stencilValues[counter] = 1.0;
          column[counter].i = i + 1;
          column[counter].j = j;
          counter++;
          counter_fluid++;
        } else {
          stencilValues[counter] = 0.0;
          column[counter].i = i + 1;
          column[counter].j = j;
          counter++;
        }
        if ((obstacle & OBSTACLE_BOTTOM) == 0) {
          stencilValues[counter] = 1.0;
          column[counter].i = i;
          column[counter].j = j - 1;
          counter++;
          counter_fluid++;
        } else {
          stencilValues[counter] = 0.0;
          column[counter].i = i;
          column[counter].j = j - 1;
          counter++;
        }
        if ((obstacle & OBSTACLE_TOP) == 0) {
          stencilValues[counter] = 1.0;
          column[counter].i = i;
          column[counter].j = j + 1;
          counter++;
          counter_fluid++;
        } else {
          stencilValues[counter] = 0.0;
          column[counter].i = i;
          column[counter].j = j + 1;
          counter++;
        }

        // A column for the cell itself
        stencilValues[counter] = -(PetscScalar)(counter_fluid);
        column[counter].i = i;
        column[counter].j = j;

        // Once we identified how many fluid cells are around and set columns
        // for each, we
        // enter the row into the matrix.

        MatSetValuesStencil(A, 1, &row, 5, column, stencilValues,
                            INSERT_VALUES);
      } else {  // The remaining possibility is that the cell is obstacle
                // surrounded
                // by more obstacle cells
        // Here, we just add an equation to set the value according to the right
        // hand side
        stencilValues[1] = 0.0;  // left
        stencilValues[0] = 0.0;  // right
        stencilValues[3] = 0.0;  // top
        stencilValues[4] = 0.0;  // bottom
        stencilValues[2] = 1.0;  // center

        // Definition of positions. Order must correspond to values
        column[0].i = i + 1;
        column[0].j = j;
        column[1].i = i - 1;
        column[1].j = j;
        column[2].i = i;
        column[2].j = j;
        column[3].i = i;
        column[3].j = j + 1;
        column[4].i = i;
        column[4].j = j - 1;
        MatSetValuesStencil(A, 1, &row, 5, column, stencilValues,
                            INSERT_VALUES);
      }
    }
  }

  // Left wall
  if (context->setAsBoundary & LEFT_WALL_BIT) {
    for (j = limitsY[0]; j < limitsY[1]; j++) {
      column[0].i = 0;
      column[0].j = j;
      column[1].i = context->displacement[0];
      column[1].j = j;
      row.i = 0;
      row.j = j;
      if (parameters.walls.typeLeft ==
          DIRICHLET) {  // If Dirichlet velocity boundary conditions
        // therefore, Neumann in the pressure
        stencilValues[0] = 1;
        stencilValues[1] = -1;
      } else if (parameters.walls.typeLeft ==
                 NEUMANN) {  // Neumann velocity boundary conditions,
        stencilValues[0] = 0.5;
        stencilValues[1] = 0.5;
      }
      MatSetValuesStencil(A, 1, &row, 2, column, stencilValues, INSERT_VALUES);
    }
  }

  // Right wall
  if (context->setAsBoundary & RIGHT_WALL_BIT) {
    for (j = limitsY[0]; j < limitsY[1]; j++) {
      column[0].i = Nx - 1;
      column[0].j = j;
      column[1].i = context->displacement[1];
      column[1].j = j;
      row.i = Nx - 1;
      row.j = j;
      if (parameters.walls.typeRight == DIRICHLET) {
        stencilValues[0] = 1;
        stencilValues[1] = -1;
      } else if (parameters.walls.typeRight == NEUMANN) {
        stencilValues[0] = 0.5;
        stencilValues[1] = 0.5;
      }
      MatSetValuesStencil(A, 1, &row, 2, column, stencilValues, INSERT_VALUES);
    }
  }

  // Bottom wall
  if (context->setAsBoundary & BOTTOM_WALL_BIT) {
    for (i = limitsX[0]; i < limitsX[1]; i++) {
      column[0].i = i;
      column[0].j = 0;
      column[1].i = i;
      column[1].j = context->displacement[2];
      row.i = i;
      row.j = 0;
      if (parameters.walls.typeBottom == DIRICHLET) {
        stencilValues[0] = 1;
        stencilValues[1] = -1;
      } else if (parameters.walls.typeBottom == NEUMANN) {
        stencilValues[0] = 0.5;
        stencilValues[1] = 0.5;
      }
      MatSetValuesStencil(A, 1, &row, 2, column, stencilValues, INSERT_VALUES);
    }
  }

  // Top wall
  if (context->setAsBoundary & TOP_WALL_BIT) {
    for (i = limitsX[0]; i < limitsX[1]; i++) {
      column[0].i = i;
      column[0].j = Ny - 1;
      column[1].i = i;
      column[1].j = context->displacement[3];
      row.i = i;
      row.j = Ny - 1;
      if (parameters.walls.typeTop == DIRICHLET) {
        stencilValues[0] = 1;
        stencilValues[1] = -1;
      } else if (parameters.walls.typeTop == NEUMANN) {
        stencilValues[0] = 0.5;
        stencilValues[1] = 0.5;
      }
      MatSetValuesStencil(A, 1, &row, 2, column, stencilValues, INSERT_VALUES);
    }
  }

  MatAssemblyBegin(A, MAT_FINAL_ASSEMBLY);
  MatAssemblyEnd(A, MAT_FINAL_ASSEMBLY);

  MatNullSpace nullspace;
  MatNullSpaceCreate(PETSC_COMM_WORLD, PETSC_TRUE, 0, 0, &nullspace);
  MatSetNullSpace(A, nullspace);
  MatNullSpaceDestroy(&nullspace);

  return 0;
}

#if (PETSC_VERSION_MAJOR == 3 && PETSC_VERSION_MINOR >= 5)
PetscErrorCode computeMatrix3D(KSP ksp, Mat A, Mat pc, void *ctx) {
#else
PetscErrorCode computeMatrix3D(KSP ksp, Mat A, Mat pc,
                               MatStructure *matStructure, void *ctx) {
#endif
  PetscUserCtx *context = (PetscUserCtx *)ctx;
  Parameters &parameters = context->getParameters();

  IntScalarField &flags = context->getFlowField().getFlags();

  int *limitsX, *limitsY, *limitsZ;
  context->getLimits(&limitsX, &limitsY, &limitsZ);

  PetscScalar stencilValues[7];
  MatStencil row, column[7];

  PetscInt i, j, k, Nx, Ny, Nz;

  Nx = parameters.geometry.sizeX + 2;
  Ny = parameters.geometry.sizeY + 2;
  Nz = parameters.geometry.sizeZ + 2;

  // Loop for inner nodes
  std::cout << "Limits: " << limitsX[0] << ", " << limitsX[1] << ", "
            << limitsY[0] << ", " << limitsY[1] << ", " << limitsZ[0] << ", "
            << limitsZ[1] << std::endl;
  for (k = limitsZ[0]; k < limitsZ[1]; k++) {
    for (j = limitsY[0]; j < limitsY[1]; j++) {
      for (i = limitsX[0]; i < limitsX[1]; i++) {
        row.i = i;
        row.j = j, row.k = k;

        const int cellIndexX = i - limitsX[0] + 2;
        const int cellIndexY = j - limitsY[0] + 2;
        const int cellIndexZ = k - limitsZ[0] + 2;
        const int obstacle = flags.getValue(cellIndexX, cellIndexY, cellIndexZ);

        if ((obstacle & OBSTACLE_SELF) == 0) {  // If the cell is fluid
          const FLOAT dx_0 =
              parameters.meshsize->getDx(cellIndexX, cellIndexY, cellIndexZ);
          const FLOAT dx_M1 = parameters.meshsize->getDx(
              cellIndexX - 1, cellIndexY, cellIndexZ);
          const FLOAT dx_P1 = parameters.meshsize->getDx(
              cellIndexX + 1, cellIndexY, cellIndexZ);
          const FLOAT dy_0 =
              parameters.meshsize->getDy(cellIndexX, cellIndexY, cellIndexZ);
          const FLOAT dy_M1 = parameters.meshsize->getDy(
              cellIndexX, cellIndexY - 1, cellIndexZ);
          const FLOAT dy_P1 = parameters.meshsize->getDy(
              cellIndexX, cellIndexY + 1, cellIndexZ);
          const FLOAT dz_0 =
              parameters.meshsize->getDz(cellIndexX, cellIndexY, cellIndexZ);
          const FLOAT dz_M1 = parameters.meshsize->getDz(cellIndexX, cellIndexY,
                                                         cellIndexZ - 1);
          const FLOAT dz_P1 = parameters.meshsize->getDz(cellIndexX, cellIndexY,
                                                         cellIndexZ + 1);

          const FLOAT dx_L = 0.5 * (dx_0 + dx_M1);
          const FLOAT dx_R = 0.5 * (dx_0 + dx_P1);
          const FLOAT dx_Bo = 0.5 * (dy_0 + dy_M1);
          const FLOAT dx_T = 0.5 * (dy_0 + dy_P1);
          const FLOAT dx_F = 0.5 * (dz_0 + dz_M1);
          const FLOAT dx_B = 0.5 * (dz_0 + dz_P1);

          // Definition of values
          stencilValues[1] = 2.0 / (dx_L * (dx_L + dx_R));    // left
          stencilValues[0] = 2.0 / (dx_R * (dx_L + dx_R));    // right
          stencilValues[2] = 2.0 / (dx_T * (dx_T + dx_Bo));   // top
          stencilValues[3] = 2.0 / (dx_Bo * (dx_T + dx_Bo));  // bottom
          stencilValues[4] = 2.0 / (dx_B * (dx_B + dx_F));    // back
          stencilValues[5] = 2.0 / (dx_F * (dx_B + dx_F));    // front
          stencilValues[6] = -2.0 / (dx_R * dx_L) - 2.0 / (dx_T * dx_Bo) -
                             2.0 / (dx_F * dx_B);  // center

          // Definition of positions. Order must correspond to values
          column[0].i = i + 1;
          column[0].j = j;
          column[0].k = k;
          column[1].i = i - 1;
          column[1].j = j;
          column[1].k = k;
          column[2].i = i;
          column[2].j = j + 1;
          column[2].k = k;
          column[3].i = i;
          column[3].j = j - 1;
          column[3].k = k;
          column[4].i = i;
          column[4].j = j;
          column[4].k = k + 1;
          column[5].i = i;
          column[5].j = j;
          column[5].k = k - 1;
          column[6].i = i;
          column[6].j = j;
          column[6].k = k;

          MatSetValuesStencil(A, 1, &row, 7, column, stencilValues,
                              INSERT_VALUES);
        } else if (obstacle !=
                   127) {  // If non-fluid and still not completely surounded
          int counter = 0;
          int counter_fluid = 0;
          if ((obstacle & OBSTACLE_LEFT) ==
              0) {  // If there's fluid to the left
            stencilValues[counter] = 1;
            column[counter].i = i - 1;
            column[counter].j = j, column[counter].k = k;
            counter++;
            counter_fluid++;
          } else {
            stencilValues[counter] = 0.0;
            column[counter].i = i - 1;
            column[counter].j = j, column[counter].k = k;
            counter++;
          }
          if ((obstacle & OBSTACLE_RIGHT) == 0) {
            stencilValues[counter] = 1;
            column[counter].i = i + 1;
            column[counter].j = j, column[counter].k = k;
            counter++;
            counter_fluid++;
          } else {
            stencilValues[counter] = 0.0;
            column[counter].i = i + 1;
            column[counter].j = j, column[counter].k = k;
            counter++;
          }
          if ((obstacle & OBSTACLE_BOTTOM) == 0) {
            stencilValues[counter] = 1;
            column[counter].i = i;
            column[counter].j = j - 1, column[counter].k = k;
            counter++;
            counter_fluid++;
          } else {
            stencilValues[counter] = 0.0;
            column[counter].i = i;
            column[counter].j = j - 1, column[counter].k = k;
            counter++;
          }
          if ((obstacle & OBSTACLE_TOP) == 0) {
            stencilValues[counter] = 1;
            column[counter].i = i;
            column[counter].j = j + 1, column[counter].k = k;
            counter++;
            counter_fluid++;
          } else {
            stencilValues[counter] = 0.0;
            column[counter].i = i;
            column[counter].j = j + 1, column[counter].k = k;
            counter++;
          }
          if ((obstacle & OBSTACLE_FRONT) == 0) {
            stencilValues[counter] = 1;
            column[counter].i = i;
            column[counter].j = j, column[counter].k = k - 1;
            counter++;
            counter_fluid++;
          } else {
            stencilValues[counter] = 0.0;
            column[counter].i = i;
            column[counter].j = j, column[counter].k = k - 1;
            counter++;
          }
          if ((obstacle & OBSTACLE_BACK) == 0) {
            stencilValues[counter] = 1;
            column[counter].i = i;
            column[counter].j = j, column[counter].k = k + 1;
            counter++;
            counter_fluid++;
          } else {
            stencilValues[counter] = 0.0;
            column[counter].i = i;
            column[counter].j = j, column[counter].k = k + 1;
            counter++;
          }

          // Now set a line for the element itself
          stencilValues[counter] = (PetscScalar)(-counter_fluid);
          column[counter].i = i;
          column[counter].j = j, column[counter].k = k;

          MatSetValuesStencil(A, 1, &row, counter + 1, column, stencilValues,
                              INSERT_VALUES);
        } else {  // The remaining possibility is that the cell is obstacle
                  // surrounded
                  // by more obstacle cells
          // Here, we just add an equation to set the value according to the
          // right hand side
          stencilValues[1] = 0.0;  // left
          stencilValues[0] = 0.0;  // right
          stencilValues[3] = 0.0;  // top
          stencilValues[4] = 0.0;  // bottom
          stencilValues[5] = 0.0;  // front
          stencilValues[5] = 0.0;  // back
          stencilValues[2] = 1.0;  // center

          // Definition of positions. Order must correspond to values
          column[0].i = i + 1;
          column[0].j = j;
          column[0].k = k;
          column[1].i = i - 1;
          column[1].j = j;
          column[1].k = k;
          column[2].i = i;
          column[2].j = j;
          column[2].k = k;
          column[3].i = i;
          column[3].j = j + 1;
          column[3].k = k;
          column[4].i = i;
          column[4].j = j - 1;
          column[4].k = k;
          column[5].i = i;
          column[5].j = j;
          column[5].k = k + 1;
          column[6].i = i;
          column[6].j = j;
          column[6].k = k - 1;
          MatSetValuesStencil(A, 1, &row, 7, column, stencilValues,
                              INSERT_VALUES);
        }
      }
    }
  }

  // Left wall
  if (context->setAsBoundary & LEFT_WALL_BIT) {
    for (j = limitsY[0]; j < limitsY[1]; j++) {
      for (k = limitsZ[0]; k < limitsZ[1]; k++) {
        column[0].i = 0;
        column[0].j = j;
        column[0].k = k;
        column[1].i = context->displacement[0];
        column[1].j = j;
        column[1].k = k;
        row.i = 0;
        row.j = j;
        row.k = k;
        if (parameters.walls.typeLeft == DIRICHLET) {
          stencilValues[0] = 1;
          stencilValues[1] = -1;
        } else if (parameters.walls.typeLeft == NEUMANN) {
          stencilValues[0] = 0.5;
          stencilValues[1] = 0.5;
        }
        MatSetValuesStencil(A, 1, &row, 2, column, stencilValues,
                            INSERT_VALUES);
      }
    }
  }

  // Right wall
  if (context->setAsBoundary & RIGHT_WALL_BIT) {
    for (j = limitsY[0]; j < limitsY[1]; j++) {
      for (k = limitsZ[0]; k < limitsZ[1]; k++) {
        column[0].i = Nx - 1;
        column[0].j = j;
        column[0].k = k;
        column[1].i = context->displacement[1];
        column[1].j = j;
        column[1].k = k;
        row.i = Nx - 1;
        row.j = j;
        row.k = k;
        if (parameters.walls.typeRight == DIRICHLET) {
          stencilValues[0] = 1;
          stencilValues[1] = -1;
        } else if (parameters.walls.typeRight == NEUMANN) {
          stencilValues[0] = 0.5;
          stencilValues[1] = 0.5;
        }
        MatSetValuesStencil(A, 1, &row, 2, column, stencilValues,
                            INSERT_VALUES);
      }
    }
  }

  // Bottom wall
  if (context->setAsBoundary & BOTTOM_WALL_BIT) {
    for (i = limitsX[0]; i < limitsX[1]; i++) {
      for (k = limitsZ[0]; k < limitsZ[1]; k++) {
        column[0].i = i;
        column[0].j = 0;
        column[0].k = k;
        column[1].i = i;
        column[1].j = context->displacement[2];
        column[1].k = k;
        row.i = i;
        row.j = 0;
        row.k = k;
        if (parameters.walls.typeBottom == DIRICHLET) {
          stencilValues[0] = 1;
          stencilValues[1] = -1;
        } else if (parameters.walls.typeBottom == NEUMANN) {
          stencilValues[0] = 0.5;
          stencilValues[1] = 0.5;
        }
        MatSetValuesStencil(A, 1, &row, 2, column, stencilValues,
                            INSERT_VALUES);
      }
    }
  }

  // Top wall
  if (context->setAsBoundary & TOP_WALL_BIT) {
    for (i = limitsX[0]; i < limitsX[1]; i++) {
      for (k = limitsZ[0]; k < limitsZ[1]; k++) {
        column[0].i = i;
        column[0].j = Ny - 1;
        column[0].k = k;
        column[1].i = i;
        column[1].j = context->displacement[3];
        column[1].k = k;
        row.i = i;
        row.j = Ny - 1;
        row.k = k;
        if (parameters.walls.typeTop == DIRICHLET) {
          stencilValues[0] = 1;
          stencilValues[1] = -1;
        } else if (parameters.walls.typeTop == NEUMANN) {
          stencilValues[0] = 0.5;
          stencilValues[1] = 0.5;
        }
        MatSetValuesStencil(A, 1, &row, 2, column, stencilValues,
                            INSERT_VALUES);
      }
    }
  }

  // Front wall
  if (context->setAsBoundary & FRONT_WALL_BIT) {
    for (i = limitsX[0]; i < limitsX[1]; i++) {
      for (j = limitsY[0]; j < limitsY[1]; j++) {
        column[0].i = i;
        column[0].j = j;
        column[0].k = 0;
        column[1].i = i;
        column[1].j = j;
        column[1].k = context->displacement[4];
        row.i = i;
        row.j = j;
        row.k = 0;
        if (parameters.walls.typeFront == DIRICHLET) {
          stencilValues[0] = 1;
          stencilValues[1] = -1;
        } else if (parameters.walls.typeFront == NEUMANN) {
          stencilValues[0] = 0.5;
          stencilValues[1] = 0.5;
        }
        MatSetValuesStencil(A, 1, &row, 2, column, stencilValues,
                            INSERT_VALUES);
      }
    }
  }

  // Back wall
  if (context->setAsBoundary & BACK_WALL_BIT) {
    for (i = limitsX[0]; i < limitsX[1]; i++) {
      for (j = limitsY[0]; j < limitsY[1]; j++) {
        column[0].i = i;
        column[0].j = j;
        column[0].k = Nz - 1;
        column[1].i = i;
        column[1].j = j;
        column[1].k = context->displacement[5];
        row.i = i;
        row.j = j;
        row.k = Nz - 1;
        if (parameters.walls.typeBack == DIRICHLET) {
          stencilValues[0] = 1;
          stencilValues[1] = -1;
        } else if (parameters.walls.typeBack == NEUMANN) {
          stencilValues[0] = 0.5;
          stencilValues[1] = 0.5;
        }
        MatSetValuesStencil(A, 1, &row, 2, column, stencilValues,
                            INSERT_VALUES);
      }
    }
  }

  MatAssemblyBegin(A, MAT_FINAL_ASSEMBLY);
  MatAssemblyEnd(A, MAT_FINAL_ASSEMBLY);

  MatNullSpace nullspace;
  MatNullSpaceCreate(PETSC_COMM_WORLD, PETSC_TRUE, 0, 0, &nullspace);
  MatSetNullSpace(A, nullspace);
  MatNullSpaceDestroy(&nullspace);

  return 0;
}

PetscErrorCode computeRHS2D(KSP ksp, Vec b, void *ctx) {
  FlowField &flowField = ((PetscUserCtx *)ctx)->getFlowField();
  Parameters &parameters = ((PetscUserCtx *)ctx)->getParameters();
  PetscUserCtx *context = ((PetscUserCtx *)ctx);

  int *limitsX, *limitsY, *limitsZ;
  ((PetscUserCtx *)ctx)->getLimits(&limitsX, &limitsY, &limitsZ);

  IntScalarField &flags = context->getFlowField().getFlags();

  ScalarField &RHS = flowField.getRHS();

  PetscInt i, j;
  PetscInt Nx = parameters.geometry.sizeX + 2,
           Ny = parameters.geometry.sizeY + 2;
  PetscScalar **array;

  DM da;
  KSPGetDM(ksp, &da);

  DMDAVecGetArray(da, b, &array);

  // Iteration domains are going to be set and the values on the global boundary
  // set when
  // necessary
  // Check left wall
  if (context->setAsBoundary & LEFT_WALL_BIT) {
    if (parameters.simulation.scenario == "pressure-channel") {
      for (j = limitsY[0]; j < limitsY[1]; j++) {
        array[j][0] = RHS.getScalar(0, j);
      }
    } else {
      for (j = limitsY[0]; j < limitsY[1]; j++) {
        array[j][0] = 0;
      }
    }
  }

  // Check right wall
  if (context->setAsBoundary & RIGHT_WALL_BIT) {
    for (j = limitsY[0]; j < limitsY[1]; j++) {
      array[j][Nx - 1] = 0;
    }
  }

  // Check bottom wall
  if (context->setAsBoundary & BOTTOM_WALL_BIT) {
    for (i = limitsX[0]; i < limitsX[1]; i++) {
      array[0][i] = 0;
    }
  }

  // Check top wall
  if (context->setAsBoundary & TOP_WALL_BIT) {
    for (i = limitsX[0]; i < limitsX[1]; i++) {
      array[Ny - 1][i] = 0;
    }
  }

  // Fill the internal nodes. We already have the values
  for (j = limitsY[0]; j < limitsY[1]; j++) {
    for (i = limitsX[0]; i < limitsX[1]; i++) {
      const int obstacle =
          flags.getValue(i - limitsX[0] + 2, j - limitsY[0] + 2);
      if ((obstacle & OBSTACLE_SELF) == 0) {  // If this is a fluid cell
        array[j][i] = RHS.getScalar(i - limitsX[0] + 2, j - limitsY[0] + 2);
      } else {
        array[j][i] = 0.0;
      }
    }
  }

  DMDAVecRestoreArray(da, b, &array);
  VecAssemblyBegin(b);
  VecAssemblyEnd(b);

  return 0;
}

PetscErrorCode computeRHS3D(KSP ksp, Vec b, void *ctx) {
  FlowField &flowField = ((PetscUserCtx *)ctx)->getFlowField();
  Parameters &parameters = ((PetscUserCtx *)ctx)->getParameters();
  ScalarField &RHS = flowField.getRHS();
  PetscUserCtx *context = ((PetscUserCtx *)ctx);

  IntScalarField &flags = flowField.getFlags();

  int *limitsX, *limitsY, *limitsZ;
  ((PetscUserCtx *)ctx)->getLimits(&limitsX, &limitsY, &limitsZ);

  PetscInt i, j, k;
  PetscInt Nx = parameters.geometry.sizeX + 2,
           Ny = parameters.geometry.sizeY + 2,
           Nz = parameters.geometry.sizeZ + 2;
  PetscScalar ***array;

  DM da;
  KSPGetDM(ksp, &da);

  DMDAVecGetArray(da, b, &array);

  // Notice that we're covering the whole surface, including corners and edges
  // Also, the actual value is taking from the parameters.

  // Left wall
  if (context->setAsBoundary & LEFT_WALL_BIT) {
    if (parameters.simulation.scenario == "pressure-channel") {
      for (k = limitsZ[0]; k < limitsZ[1]; k++) {
        for (j = limitsY[0]; j < limitsY[1]; j++) {
          array[k][j][0] = RHS.getScalar(0, j, k);
        }
      }
    } else {
      for (k = limitsZ[0]; k < limitsZ[1]; k++) {
        for (j = limitsY[0]; j < limitsY[1]; j++) {
          array[k][j][0] = 0;
        }
      }
    }
  }

  // Right wall
  if (context->setAsBoundary & RIGHT_WALL_BIT) {
    for (k = limitsZ[0]; k < limitsZ[1]; k++) {
      for (j = limitsY[0]; j < limitsY[1]; j++) {
        array[k][j][Nx - 1] = 0;
      }
    }
  }

  // Bottom wall
  if (context->setAsBoundary & BOTTOM_WALL_BIT) {
    for (k = limitsZ[0]; k < limitsZ[1]; k++) {
      for (i = limitsX[0]; i < limitsX[1]; i++) {
        array[k][0][i] = 0;
      }
    }
  }

  // Top wall
  if (context->setAsBoundary & TOP_WALL_BIT) {
    for (k = limitsZ[0]; k < limitsZ[1]; k++) {
      for (i = limitsX[0]; i < limitsX[1]; i++) {
        array[k][Ny - 1][i] = 0;
      }
    }
  }

  // Front wall
  if (context->setAsBoundary & FRONT_WALL_BIT) {
    for (j = limitsY[0]; j < limitsY[1]; j++) {
      for (i = limitsX[0]; i < limitsX[1]; i++) {
        array[0][j][i] = 0;
      }
    }
  }

  // Back wall
  if (context->setAsBoundary & BACK_WALL_BIT) {
    for (j = limitsY[0]; j < limitsY[1]; j++) {
      for (i = limitsX[0]; i < limitsX[1]; i++) {
        array[Nz - 1][j][i] = 0;
      }
    }
  }

  // Fill the internal nodes. We already have the values
  for (k = limitsZ[0]; k < limitsZ[1]; k++) {
    for (j = limitsY[0]; j < limitsY[1]; j++) {
      for (i = limitsX[0]; i < limitsX[1]; i++) {
        const int obstacle = flags.getValue(
            i - limitsX[0] + 2, j - limitsY[0] + 2, k - limitsZ[0] + 2);
        if ((obstacle & OBSTACLE_SELF) == 0) {  // If this is a fluid cell
          array[k][j][i] = RHS.getScalar(i - limitsX[0] + 2, j - limitsY[0] + 2,
                                         k - limitsZ[0] + 2);
        } else {
          array[k][j][i] = 0.0;
        }
      }
    }
  }

  DMDAVecRestoreArray(da, b, &array);
  VecAssemblyBegin(b);
  VecAssemblyEnd(b);

  return 0;
}

const DM &PetscSolver::getGrid() const { return _da; }

void PetscUserCtx::setLimits(int *limitsX, int *limitsY, int *limitsZ) {
  _limitsX = limitsX;
  _limitsY = limitsY;
  _limitsZ = limitsZ;
}

void PetscUserCtx::getLimits(int **limitsX, int **limitsY, int **limitsZ) {
  *limitsX = _limitsX;
  *limitsY = _limitsY;
  *limitsZ = _limitsZ;
}

void PetscUserCtx::setRank(int rank) { _rank = rank; }

int PetscUserCtx::getRank() const { return _rank; }

void PetscSolver::reInitMatrix() {
  std::cout << "Reinit the matrix" << std::endl;
  if (_parameters.geometry.dim == 2)
    KSPSetComputeOperators(_ksp, computeMatrix2D, &_ctx);
  else
    KSPSetComputeOperators(_ksp, computeMatrix3D, &_ctx);
}
}
