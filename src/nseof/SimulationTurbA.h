#ifndef _SIMULATION_TURB_A_H_
#define _SIMULATION_TURB_A_H_

#include <petscksp.h>
#include <float.h>
#include "FlowField.h"
#include "FlowFieldTurbA.h"
#include "stencils/FGHStencilTurb.h"
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
#include "stencils/NutStencil.h"
#include "stencils/HStencil.h"
#include "stencils/MinimumNutStencil.h"
#include "GlobalBoundaryFactory.h"
#include "Iterators.h"
#include "Definitions.h"
#include "Simulation.h"

#include "LinearSolver.h"
#include "solvers/SORSolver.h"
#include "solvers/PetscSolver.h"
#include "parallelManagers/MPICommunicator.h"

namespace nseof {

class SimulationTurbA : public FlowFieldSimulation<FlowFieldTurbA> {
 protected:
  MaxUStencil _maxUStencil;
  FieldIterator<FlowField> _maxUFieldIterator;
  GlobalBoundaryIterator<FlowField> _maxUBoundaryIterator;

  // Set up the boundary conditions
  GlobalBoundaryFactory _globalBoundaryFactory;
  GlobalBoundaryIterator<FlowField> _wallVelocityIterator;
  GlobalBoundaryIterator<FlowField> _wallFGHIterator;

  FGHStencilTurb _fghStencil;
  FieldIterator<FlowFieldTurbA> _fghIterator;

  RHSStencil _rhsStencil;
  FieldIterator<FlowField> _rhsIterator;

  VelocityStencil _velocityStencil;
  ObstacleStencil _obstacleStencil;
  FieldIterator<FlowField> _velocityIterator;
  FieldIterator<FlowField> _obstacleIterator;

  PetscSolver _solver;

  NutStencil _nutst;
  FieldIterator<FlowFieldTurbA> _nutit;
  HStencil _hst;
  FieldIterator<FlowFieldTurbA> _hit;

  MinimumNutStencil _minnutst;
  FieldIterator<FlowFieldTurbA> _minnutit;

  MPICommunicator<FLOAT, FlowFieldTurbA> nutComm{
      *this->_flowField, this->_parameters,
      [](FlowFieldTurbA &flowField, int i, int j, int k, FLOAT &p) {
        p = flowField.getNu(i, j, k);
      },
      [](FlowFieldTurbA &flowField, int i, int j, int k, FLOAT &p) {
        flowField.getNu(i, j, k) = p;
      },
      2};

 public:
  SimulationTurbA(Parameters &parameters)
      : FlowFieldSimulation(parameters, new FlowFieldTurbA(parameters)),
        _maxUStencil(parameters),
        _maxUFieldIterator(*_flowField, parameters, _maxUStencil),
        _maxUBoundaryIterator(*_flowField, parameters, _maxUStencil),
        _globalBoundaryFactory(parameters),
        _wallVelocityIterator(
            _globalBoundaryFactory.getGlobalBoundaryVelocityIterator(
                *_flowField)),
        _wallFGHIterator(
            _globalBoundaryFactory.getGlobalBoundaryFGHIterator(*_flowField)),
        _fghStencil(parameters),
        _fghIterator(*_flowField, parameters, _fghStencil),
        _rhsStencil(parameters),
        _rhsIterator(*_flowField, parameters, _rhsStencil),
        _velocityStencil(parameters),
        _obstacleStencil(parameters),
        _velocityIterator(*_flowField, parameters, _velocityStencil),
        _obstacleIterator(*_flowField, parameters, _obstacleStencil),
        _solver(*_flowField, parameters),
        _nutst(parameters),
        _nutit(*_flowField, _parameters, _nutst, 1, 0),
        _hst(parameters),
        _hit(*_flowField, _parameters, _hst, 0, 0),
        _minnutst(parameters),
        _minnutit(*_flowField, _parameters, _minnutst, 1, 0) {
    // distance to the next wall
    this->scalarStencils.push_back(CellDataStencil<double, FlowFieldTurbA>(
        this->_parameters, "h",
        [](FlowFieldTurbA &f, int i, int j) { return f.getH(i, j); },
        [](FlowFieldTurbA &f, int i, int j, int k) {
          return f.getH(i, j, k);
        }));

    // vortex viscosity
    this->scalarStencils.push_back(CellDataStencil<double, FlowFieldTurbA>(
        this->_parameters, "nu",
        [](FlowFieldTurbA &f, int i, int j) { return f.getNu(i, j); },
        [](FlowFieldTurbA &f, int i, int j, int k) {
          return f.getNu(i, j, k);
        }));

    // mixing length
    this->scalarStencils.push_back(CellDataStencil<double, FlowFieldTurbA>(
        this->_parameters, "lm",
        [](FlowFieldTurbA &f, int i, int j) { return f.getLm(i, j); },
        [](FlowFieldTurbA &f, int i, int j, int k) {
          return f.getLm(i, j, k);
        }));

    // velocity fluctuation
    this->scalarStencils.push_back(CellDataStencil<double, FlowFieldTurbA>(
        this->_parameters, "u",
        [](FlowFieldTurbA &f, int i, int j) { return f.getU(i, j); },
        [](FlowFieldTurbA &f, int i, int j, int k) {
          return f.getU(i, j, k);
        }));

    // actual pressure (without tke)
    this->scalarStencils.push_back(CellDataStencil<double, FlowFieldTurbA>(
        this->_parameters, "pressureactual",
        [](FlowFieldTurbA &f, int i, int j) {
          return f.getPressure().getScalar(i, j) - f.getU(i, j) * f.getU(i, j);
        },
        [](FlowFieldTurbA &f, int i, int j, int k) {
          return f.getPressure().getScalar(i, j, k) -
                 f.getU(i, j, k) * f.getU(i, j, k);
        }));
  }

  virtual ~SimulationTurbA() {}

  /** initialises the flow field according to the scenario */
  virtual void initializeFlowField() {
    if (_parameters.simulation.scenario == "taylor-green") {
      // currently, a particular initialization is only requrid for the
      // taylor-green vortex
      InitTaylorGreenFlowFieldStencil stencil(_parameters);
      FieldIterator<FlowField> iterator(*_flowField, _parameters, stencil);
      iterator.iterate();
    } else if (_parameters.simulation.scenario == "channel") {
      BFStepInitStencil stencil(_parameters);
      FieldIterator<FlowField> iterator(*_flowField, _parameters, stencil, 0,
                                        1);
      iterator.iterate();
      _wallVelocityIterator.iterate();
    } else if (_parameters.simulation.scenario == "pressure-channel") {
      // set pressure boundaries here for left wall
      const FLOAT value = _parameters.walls.scalarLeft;
      ScalarField &rhs = _flowField->getRHS();

      if (_parameters.geometry.dim == 2) {
        const int sizey = _flowField->getNy();
        for (int i = 0; i < sizey + 3; i++) {
          rhs.getScalar(0, i) = value;
        }
      } else {
        const int sizey = _flowField->getNy();
        const int sizez = _flowField->getNz();
        for (int i = 0; i < sizey + 3; i++)
          for (int j = 0; j < sizez + 3; j++) rhs.getScalar(0, i, j) = value;
      }

      // do same procedure for domain flagging as for regular channel
      BFStepInitStencil stencil(_parameters);
      FieldIterator<FlowField> iterator(*_flowField, _parameters, stencil, 0,
                                        1);
      iterator.iterate();
    }
    _solver.reInitMatrix();

    _hit.iterate();
  }

  virtual void solveTimestep() {
    MultiTimer *timer = MultiTimer::get();

    // determine and set max. timestep which is allowed in this simulation
    setTimeStep();

    timer->start("fgh");

    // compute fgh (turbulent)
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
    pressureComm.communicate(*this->_flowField);

    timer->stop("pressure-communication");
    timer->stop("communication");

    // compute velocity
    _velocityIterator.iterate();
    // set obstacle boundaries
    _obstacleIterator.iterate();

    timer->start("communication");
    timer->start("velocity-communication");

    // TODO WS2: communicate velocity values
    velocityComm.communicate(*this->_flowField);

    timer->stop("velocity-communication");
    timer->stop("communication");

    // Iterate for velocities on the boundary
    _wallVelocityIterator.iterate();

    // compute vortex viscosity
    _nutit.iterate();

    // communicate vortex viscosity
    nutComm.communicate(*this->_flowField);
  }

 protected:
  /** sets the time step*/
  virtual void setTimeStep() {
    _minnutst.reset();
    _minnutit.iterate();

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

    // determin minimal timestep on domain with formula:
    // dt = 1/(2 * (nu+nut))/(dx_min^-2+dy_min^-2+dz_min^-2))
    localMin = std::min(_parameters.timestep.dt,
                        std::min(std::min(1 / (1 / _parameters.flow.Re +
                                               _minnutst.getMinimum()) /
                                              (2 * factor),
                                          1.0 / _maxUStencil.getMaxValues()[0]),
                                 1.0 / _maxUStencil.getMaxValues()[1]));

    globalMin = MY_FLOAT_MAX;
    MPI_Allreduce(&localMin, &globalMin, 1, MY_MPI_FLOAT, MPI_MIN,
                  PETSC_COMM_WORLD);

    _parameters.timestep.dt = globalMin;
    _parameters.timestep.dt *= _parameters.timestep.tau;
  }
};
}

#endif  // _SIMULATION_TURB_A_H_
