#ifndef _SIMULATION_LAMINAR_H_
#define _SIMULATION_LAMINAR_H_

#include <petscksp.h>
#include <float.h>
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
#include "GlobalBoundaryFactory.h"
#include "Iterators.h"
#include "Definitions.h"
#include "FlowFieldSimulation.h"

#include "LinearSolver.h"
#include "solvers/SORSolver.h"
#include "solvers/PetscSolver.h"

class SimulationLaminar : public FlowFieldSimulation<FlowField> {
 protected:
  MaxUStencil _maxUStencil;
  FieldIterator<FlowField> _maxUFieldIterator;
  GlobalBoundaryIterator<FlowField> _maxUBoundaryIterator;

  // Set up the boundary conditions
  GlobalBoundaryFactory _globalBoundaryFactory;
  GlobalBoundaryIterator<FlowField> _wallVelocityIterator;
  GlobalBoundaryIterator<FlowField> _wallFGHIterator;

  FGHStencil _fghStencil;
  FieldIterator<FlowField> _fghIterator;

  RHSStencil _rhsStencil;
  FieldIterator<FlowField> _rhsIterator;

  VelocityStencil _velocityStencil;
  ObstacleStencil _obstacleStencil;
  FieldIterator<FlowField> _velocityIterator;
  FieldIterator<FlowField> _obstacleIterator;

  PetscSolver _solver;

 public:
  SimulationLaminar(Parameters &parameters, FlowField &flowField)
      : FlowFieldSimulation(parameters, flowField),
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

  virtual ~SimulationLaminar() {}

  /** initialises the flow field according to the scenario */
  virtual void initializeFlowField() {
    if (_parameters.simulation.scenario == "taylor-green") {
      // currently, a particular initialization is only requrid for the
      // taylor-green vortex
      InitTaylorGreenFlowFieldStencil stencil(_parameters);
      FieldIterator<FlowField> iterator(_flowField, _parameters, stencil);
      iterator.iterate();
    } else if (_parameters.simulation.scenario == "channel") {
      BFStepInitStencil stencil(_parameters);
      FieldIterator<FlowField> iterator(_flowField, _parameters, stencil, 0, 1);
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
      FieldIterator<FlowField> iterator(_flowField, _parameters, stencil, 0, 1);
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

    localMin = 1e-4;

    // Here, we select the type of operation before compiling. This allows to
    // use the correct
    // data type for MPI. Not a concern for small simulations, but useful if
    // using heterogeneous
    // machines.

    globalMin = MY_FLOAT_MAX;
    MPI_Allreduce(&localMin, &globalMin, 1, MY_MPI_FLOAT, MPI_MIN,
                  PETSC_COMM_WORLD);

    _parameters.timestep.dt = globalMin;
    _parameters.timestep.dt = 1e-3;
    _parameters.timestep.dt *= _parameters.timestep.tau;
  }
};

#endif  // _SIMULATION_LAMINAR_H_
