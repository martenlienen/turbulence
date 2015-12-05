#ifndef _FLOW_FIELD_SIMULATION_H_
#define _FLOW_FIELD_SIMULATION_H_

#include <algorithm>
#include <array>
#include <memory>
#include <utility>
#include <vector>

#include <petscksp.h>
#include <float.h>
#include <vtk/File.h>
#include <vtk/Dataset.h>
#include <vtk/CellData.h>

#include "Simulation.h"
#include "FlowField.h"

#include "stencils/FGHStencil.h"
#include "stencils/MovingWallStencils.h"
#include "stencils/RHSStencil.h"
#include "stencils/VelocityStencil.h"
#include "stencils/ObstacleStencil.h"
#include "stencils/MaxUStencil.h"
#include "stencils/PeriodicBoundaryStencils.h"
#include "stencils/BFStepInitStencil.h"
#include "stencils/NeumannBoundaryStencils.h"
#include "stencils/BFInputStencils.h"
#include "stencils/InitTaylorGreenFlowFieldStencil.h"
#include "stencils/CellDataStencil.h"

#include "GlobalBoundaryFactory.h"
#include "Iterators.h"
#include "Definitions.h"

#include "LinearSolver.h"
#include "solvers/SORSolver.h"
#include "solvers/PetscSolver.h"

#include "parallelManagers/MPICommunicator.h"

#include "MultiTimer.h"

#define GHOST_OFFSET 2

template <typename FF>
class FlowFieldSimulation : public Simulation {
 protected:
  FF &_flowField;

  MaxUStencil _maxUStencil;
  FieldIterator<FF> _maxUFieldIterator;
  GlobalBoundaryIterator<FF> _maxUBoundaryIterator;

  // Set up the boundary conditions
  GlobalBoundaryFactory _globalBoundaryFactory;
  GlobalBoundaryIterator<FF> _wallVelocityIterator;
  GlobalBoundaryIterator<FF> _wallFGHIterator;

  FGHStencil _fghStencil;
  FieldIterator<FF> _fghIterator;

  RHSStencil _rhsStencil;
  FieldIterator<FF> _rhsIterator;

  VelocityStencil _velocityStencil;
  ObstacleStencil _obstacleStencil;
  FieldIterator<FF> _velocityIterator;
  FieldIterator<FF> _obstacleIterator;

  MPICommunicator<FLOAT, FF> pressureComm{
      this->_flowField, this->_parameters,
      [](FF &flowField, int i, int j, int k, FLOAT &p) {
        p = flowField.getPressure().getScalar(i, j, k);
      },
      [](FF &flowField, int i, int j, int k, FLOAT &p) {
        flowField.getPressure().getScalar(i, j, k) = p;
      }};

  MPICommunicator<std::array<FLOAT, 3>, FF> velocityComm{
      this->_flowField, this->_parameters,
      [](FF &flowField, int i, int j, int k, std::array<FLOAT, 3> &v) {
        std::copy_n(flowField.getVelocity().getVector(i, j, k), 3, v.data());
      },
      [](FF &flowField, int i, int j, int k, std::array<FLOAT, 3> &v) {
        std::copy_n(v.data(), 3, flowField.getVelocity().getVector(i, j, k));
      },
      2};

  std::vector<CellDataStencil<double, FF>> scalarStencils{
      CellDataStencil<double, FF>(this->_parameters, "pressure",
                                  [](FF &f, int i, int j) {
                                    FLOAT p;
                                    FLOAT v[3];
                                    f.getPressureAndVelocity(p, v, i, j);
                                    return p;
                                  },
                                  [](FF &f, int i, int j, int k) {
                                    FLOAT p;
                                    FLOAT v[3];
                                    f.getPressureAndVelocity(p, v, i, j, k);
                                    return p;
                                  })};

  std::vector<CellDataStencil<std::vector<double>, FF>> vectorStencils{
      CellDataStencil<std::vector<double>, FF>(
          this->_parameters, "velocity",
          [](FF &f, int i, int j) {
            FLOAT p;
            FLOAT v[3];
            f.getPressureAndVelocity(p, v, i, j);
            return std::vector<double>{v[0], v[1], v[2]};
          },
          [](FF &f, int i, int j, int k) {
            FLOAT p;
            FLOAT v[3];
            f.getPressureAndVelocity(p, v, i, j, k);
            return std::vector<double>{v[0], v[1], v[2]};
          })};

  PetscSolver _solver;

 public:
  FlowFieldSimulation(Parameters &parameters, FF &flowField)
      : Simulation(parameters),
        _flowField(flowField),
        _maxUStencil(parameters),
        _maxUFieldIterator(_flowField, parameters, _maxUStencil),
        _maxUBoundaryIterator(_flowField, parameters, _maxUStencil),
        _globalBoundaryFactory(parameters),
        _wallVelocityIterator(
            _globalBoundaryFactory.getGlobalBoundaryVelocityIterator(
                _flowField)),
        _wallFGHIterator(
            _globalBoundaryFactory.getGlobalBoundaryFGHIterator(_flowField)),
        _fghStencil(parameters),
        _fghIterator(_flowField, parameters, _fghStencil),
        _rhsStencil(parameters),
        _rhsIterator(_flowField, parameters, _rhsStencil),
        _velocityStencil(parameters),
        _obstacleStencil(parameters),
        _velocityIterator(_flowField, parameters, _velocityStencil),
        _obstacleIterator(_flowField, parameters, _obstacleStencil),
        _solver(_flowField, parameters) {}

  virtual ~FlowFieldSimulation() {}

  /** initialises the flow field according to the scenario */
  virtual void initializeFlowField() {
    if (_parameters.simulation.scenario == "taylor-green") {
      // currently, a particular initialization is only requrid for the
      // taylor-green vortex
      InitTaylorGreenFlowFieldStencil stencil(_parameters);
      FieldIterator<FF> iterator(_flowField, _parameters, stencil);
      iterator.iterate();
    } else if (_parameters.simulation.scenario == "channel") {
      BFStepInitStencil stencil(_parameters);
      FieldIterator<FF> iterator(_flowField, _parameters, stencil, 0, 1);
      iterator.iterate();
      _wallVelocityIterator.iterate();
    } else if (_parameters.simulation.scenario == "pressure-channel") {
      // set pressure boundaries here for left wall
      const FLOAT value = _parameters.walls.scalarLeft;
      ScalarField &rhs = _flowField.getRHS();

      if (_parameters.geometry.dim == 2) {
        const int sizey = _flowField.getNy();
        for (int i = 0; i < sizey + 3; i++) {
          rhs.getScalar(0, i) = value;
        }
      } else {
        const int sizey = _flowField.getNy();
        const int sizez = _flowField.getNz();
        for (int i = 0; i < sizey + 3; i++)
          for (int j = 0; j < sizez + 3; j++) rhs.getScalar(0, i, j) = value;
      }

      // do same procedure for domain flagging as for regular channel
      BFStepInitStencil stencil(_parameters);
      FieldIterator<FF> iterator(_flowField, _parameters, stencil, 0, 1);
      iterator.iterate();
    }
    _solver.reInitMatrix();
  }

  virtual void solveTimestep() {
    MultiTimer *timer = MultiTimer::get();

    // determine and set max. timestep which is allowed in this simulation
    setTimeStep();

    timer->start("fgh");

    // compute fgh
    _fghIterator.iterate();

    timer->stop("fgh");

    // set global boundary values
    _wallFGHIterator.iterate();
    // compute the right hand side
    _rhsIterator.iterate();

    timer->start("poisson");

    // solve for pressure
    _solver.solve();

    timer->stop("poisson");

    timer->start("communication");
    timer->start("pressure-communication");

    // TODO WS2: communicate pressure values
    pressureComm.communicate(this->_flowField);

    timer->stop("pressure-communication");
    timer->stop("communication");

    // compute velocity
    _velocityIterator.iterate();
    // set obstacle boundaries
    _obstacleIterator.iterate();

    timer->start("communication");
    timer->start("velocity-communication");

    // TODO WS2: communicate velocity values
    velocityComm.communicate(this->_flowField);

    timer->stop("velocity-communication");
    timer->stop("communication");

    // Iterate for velocities on the boundary
    _wallVelocityIterator.iterate();
  }

  /** TODO WS1: plots the flow field. */
  virtual void plotVTK(int rank, int timeStep) {
    if (!this->_parameters.vtk.enabled) {
      return;
    }

    // TODO WS1: create VTKStencil and respective iterator; iterate stencil
    //           over _flowField and write flow field information to vtk file

    vtk::Dataset dataset = this->datasetFromMesh();
    std::vector<std::unique_ptr<vtk::CellData>> cellData;
    cellData.reserve(this->scalarStencils.size() + this->vectorStencils.size());

    for (auto &stencil : this->scalarStencils) {
      stencil.reset();
      FieldIterator<FF> iterator(this->_flowField, this->_parameters, stencil,
                                 1, 0);
      iterator.iterate();
      cellData.push_back(stencil.get());
    }

    for (auto &stencil : this->vectorStencils) {
      stencil.reset();
      FieldIterator<FF> iterator(this->_flowField, this->_parameters, stencil,
                                 1, 0);
      iterator.iterate();
      cellData.push_back(stencil.get());
    }

    vtk::File file(dataset, std::move(cellData));
    file.write(this->_parameters.vtk.prefix, rank, timeStep);
  }

 protected:
  /** sets the time step*/
  virtual void setTimeStep() {
    FLOAT localMin, globalMin;
    assertion(_parameters.geometry.dim == 2 || _parameters.geometry.dim == 3);
    FLOAT factor = 1.0 / (_parameters.meshsize->getDxMin() *
                          _parameters.meshsize->getDxMin()) +
                   1.0 / (_parameters.meshsize->getDyMin() *
                          _parameters.meshsize->getDyMin());

    // determine maximum velocity
    _maxUStencil.reset();
    _maxUFieldIterator.iterate();
    _maxUBoundaryIterator.iterate();
    if (_parameters.geometry.dim == 3) {
      factor += 1.0 / (_parameters.meshsize->getDzMin() *
                       _parameters.meshsize->getDzMin());
      _parameters.timestep.dt = 1.0 / _maxUStencil.getMaxValues()[2];
    } else {
      _parameters.timestep.dt = 1.0 / _maxUStencil.getMaxValues()[0];
    }

    localMin = std::min(_parameters.timestep.dt,
                        std::min(std::min(_parameters.flow.Re / (2 * factor),
                                          1.0 / _maxUStencil.getMaxValues()[0]),
                                 1.0 / _maxUStencil.getMaxValues()[1]));

    // Here, we select the type of operation before compiling. This allows to
    // use the correct
    // data type for MPI. Not a concern for small simulations, but useful if
    // using heterogeneous
    // machines.

    globalMin = MY_FLOAT_MAX;
    MPI_Allreduce(&localMin, &globalMin, 1, MY_MPI_FLOAT, MPI_MIN,
                  PETSC_COMM_WORLD);

    _parameters.timestep.dt = globalMin;
    _parameters.timestep.dt *= _parameters.timestep.tau;
  }

 private:
  vtk::Dataset datasetFromMesh() {
    Parameters &p = this->_parameters;
    ParallelParameters &pp = p.parallel;
    GeometricParameters &gp = p.geometry;
    Meshsize *mesh = p.meshsize;

    int cellsX = pp.localSize[0];
    int cellsY = pp.localSize[1];
    int cellsZ = gp.dim == 3 ? pp.localSize[2] : 1;
    int pointsX = cellsX + 1;
    int pointsY = cellsY + 1;
    int pointsZ = gp.dim == 3 ? cellsZ + 1 : 1;
    int numPoints = pointsX * pointsY * pointsZ;

    std::vector<vtk::Point> points(numPoints);

    for (int k = GHOST_OFFSET, p = 0; k < pointsZ + GHOST_OFFSET; k++) {
      for (int j = GHOST_OFFSET; j < pointsY + GHOST_OFFSET; j++) {
        for (int i = GHOST_OFFSET; i < pointsX + GHOST_OFFSET; i++, p++) {
          points[p].x = mesh->getPosX(i, j, k);
          points[p].y = mesh->getPosY(i, j, k);
          points[p].z = mesh->getPosZ(i, j, k);
        }
      }
    }

    return vtk::Dataset(vtk::Point(cellsX, cellsY, cellsZ),
                        vtk::Point(pointsX, pointsY, pointsZ), points);
  }
};

#endif  // _FLOW_FIELD_SIMULATION_H_
