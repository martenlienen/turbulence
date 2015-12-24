#include <stdlib.h>
#include <stdio.h>
#include <iostream>
#include <fstream>
#include <sstream>
#include <iomanip>

#include "Configuration.h"
#include "Simulation.h"
#include "SimulationLaminar.h"
#include "SimulationTurbA.h"
#include "parallelManagers/PetscParallelConfiguration.h"
#include "MeshsizeFactory.h"
#include "FlowFieldTurbA.h"
#include "FlowFieldLaminar.h"
#include "MultiTimer.h"

int main(int argc, char *argv[]) {
  auto timer = nseof::MultiTimer::get();
  timer->start("total");
  timer->start("initialization");

  // Parallelization related. Initialize and identify
  // ---------------------------------------------------
  int rank;   // This processor's identifier
  int nproc;  // Number of processors in the group
  PetscInitialize(&argc, &argv, "petsc_commandline_arg", PETSC_NULL);
  MPI_Comm_size(PETSC_COMM_WORLD, &nproc);
  MPI_Comm_rank(PETSC_COMM_WORLD, &rank);
  std::cout << "Rank: " << rank << ", Nproc: " << nproc << std::endl;
  //----------------------------------------------------
  // read configuration and store information in parameters object
  nseof::Configuration configuration(argv[1]);
  nseof::Parameters parameters;
  configuration.loadParameters(parameters);
  nseof::PetscParallelConfiguration parallelConfiguration(parameters);
  nseof::MeshsizeFactory::getInstance().initMeshsize(parameters);
  nseof::FlowField *flowField = NULL;
  nseof::Simulation *simulation = NULL;

  if (parameters.geometry.dim == 2) {
    parameters.geometry.sizeZ = 1;
    parameters.parallel.localSize[2] = 1;
  }

#ifdef DEBUG
  std::cout << "Processor " << parameters.parallel.rank << " with index ";
  std::cout << parameters.parallel.indices[0] << ",";
  std::cout << parameters.parallel.indices[1] << ",";
  std::cout << parameters.parallel.indices[2];
  std::cout << " is computing the size of its subdomain and obtains ";
  std::cout << parameters.parallel.localSize[0] << ", ";
  std::cout << parameters.parallel.localSize[1] << " and ";
  std::cout << parameters.parallel.localSize[2];
  std::cout << ". Left neighbour: " << parameters.parallel.leftNb;
  std::cout << ", right neighbour: " << parameters.parallel.rightNb;
  std::cout << std::endl;
  std::cout << "Min. meshsizes: " << parameters.meshsize->getDxMin() << ", "
            << parameters.meshsize->getDyMin() << ", "
            << parameters.meshsize->getDzMin() << std::endl;
#endif

  // initialise simulation
  if (parameters.simulation.type == "turbulence") {
    if (rank == 0) {
      std::cout << "Start turbulent simulation in " << parameters.geometry.dim
                << "D" << std::endl;
    }

    // create algebraic turbulent flow field
    auto flowFieldT = new nseof::FlowFieldTurbA(parameters);
    if (flowFieldT == NULL) {
      handleError(1, "flowFieldT==NULL!");
    }

    // create algebraic turbulent simulation
    simulation = new nseof::SimulationTurbA(parameters, *flowFieldT);

    flowField = flowFieldT;
  } else if (parameters.simulation.type == "dns") {
    if (rank == 0) {
      std::cout << "Start DNS simulation in " << parameters.geometry.dim << "D"
                << std::endl;
    }

    flowField = new nseof::FlowFieldLaminar(parameters);
    if (flowField == NULL) {
      handleError(1, "flowField==NULL!");
    }

    simulation = new nseof::SimulationLaminar(parameters, *flowField);
  } else {
    handleError(
        1, "Unknown simulation type! Currently supported: dns, turbulence");
  }
  // call initialization of simulation (initialize flow field)
  if (simulation == NULL) {
    handleError(1, "simulation==NULL!");
  }
  simulation->initializeFlowField();
  // flowField->getFlags().show();

  nseof::FLOAT time = 0.0;
  nseof::FLOAT timeVTK = parameters.vtk.interval;
  int timeSteps = 0;

  // TODO WS1: plot initial state
  simulation->plotVTK(rank, 0);

  timer->stop("initialization");

  // time loop
  while (time < parameters.simulation.finalTime) {
    simulation->solveTimestep();

    time += parameters.timestep.dt;

    // std-out: terminal info
    if (rank == 0) {
      std::cout << "Current time: " << time
                << "\ttimestep: " << parameters.timestep.dt << std::endl;
    }

    timeSteps++;

    // TODO WS1: trigger VTK output
    if (time >= timeVTK) {
      simulation->plotVTK(rank, timeSteps);
      timeVTK += parameters.vtk.interval;
    }
  }

  // TODO WS1: plot final output
  simulation->plotVTK(rank, timeSteps);

  delete simulation;
  simulation = NULL;
  delete flowField;
  flowField = NULL;

  PetscFinalize();

  timer->stop("total");

  if (parameters.timing.enabled) {
    std::ostringstream path;
    path << parameters.timing.prefix << "rank-" << rank;

    std::fstream file;
    file.open(path.str().c_str(), std::ios::out);

    if (file.fail()) {
      std::cout << "Could not open " << std::endl;
      return 1;
    }

    file << timer->toString();
    file << "timesteps " << timeSteps << std::endl;
  }
}
