#ifndef _NSEOF_FLOWMODELS_ALGEBRAIC_SIMULATION_H_
#define _NSEOF_FLOWMODELS_ALGEBRAIC_SIMULATION_H_

#include <algorithm>

#include <petscksp.h>
#include <float.h>

#include "../../FlowField.h"
#include "../../GlobalBoundaryFactory.h"
#include "../../Iterators.h"
#include "../../Definitions.h"
#include "../../FlowFieldSimulation.h"
#include "../../stencils/RHSStencil.h"
#include "../../stencils/VelocityStencil.h"
#include "../../stencils/ObstacleStencil.h"
#include "../../stencils/MaxUStencil.h"
#include "../../stencils/BFStepInitStencil.h"
#include "../../stencils/InitTaylorGreenFlowFieldStencil.h"
#include "../../solvers/PetscSolver.h"
#include "../../parallelManagers/MPICommunicator.h"

#include "../turbulent/FlowField.h"
#include "../turbulent/HStencil.h"
#include "../turbulent/MaximumNutStencil.h"
#include "../turbulent/FGHStencil.h"

#include "FlowField.h"
#include "NutStencil.h"

namespace nseof {

namespace flowmodels {

namespace algebraic {

class Simulation : public FlowFieldSimulation<FlowField> {
 protected:
  MaxUStencil _maxUStencil;
  FieldIterator<nseof::FlowField> _maxUFieldIterator;
  GlobalBoundaryIterator<nseof::FlowField> _maxUBoundaryIterator;

  // Set up the boundary conditions
  GlobalBoundaryFactory _globalBoundaryFactory;
  GlobalBoundaryIterator<nseof::FlowField> _wallVelocityIterator;
  GlobalBoundaryIterator<nseof::FlowField> _wallFGHIterator;

  nseof::flowmodels::turbulent::FGHStencil _fghStencil;
  FieldIterator<nseof::flowmodels::turbulent::FlowField> _fghIterator;

  RHSStencil _rhsStencil;
  FieldIterator<nseof::FlowField> _rhsIterator;

  VelocityStencil _velocityStencil;
  ObstacleStencil _obstacleStencil;
  FieldIterator<nseof::FlowField> _velocityIterator;
  FieldIterator<nseof::FlowField> _obstacleIterator;

  std::vector<std::string> _mpiiwvector2D;
  std::vector<std::string> _mpiiwvector3D;
  MPIIteratorWrite<FlowField, double> _mpiiw;
  MPIIteratorRead<FlowField, double> _mpiir;

  NutStencil _nutst;
  FieldIterator<FlowField> _nutit;
  nseof::flowmodels::turbulent::HStencil _hst;
  FieldIterator<nseof::flowmodels::turbulent::FlowField> _hit;

  nseof::flowmodels::turbulent::MaximumNutStencil _maxmnutst;
  FieldIterator<nseof::flowmodels::turbulent::FlowField> _maxmnutit;

  MPICommunicator<FLOAT, FlowField> nutComm{
      *this->_flowField, this->_parameters,
      [](FlowField &flowField, int i, int j, int k, FLOAT &p) {
        p = flowField.getNu(i, j, k);
      },
      [](FlowField &flowField, int i, int j, int k, FLOAT &p) {
        flowField.getNu(i, j, k) = p;
      },
      2};

 public:
  Simulation(Parameters &parameters, nseof::geometry::GeometryManager &gm)
      : FlowFieldSimulation(parameters, new FlowField(parameters), gm),
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
        _mpiiwvector2D{"p", "u", "v", "nut", "P"},
        _mpiiwvector3D{"p", "u", "v", "w", "nut", "P"},
        _mpiiw(
            *_flowField, parameters, _mpiiwvector2D, _mpiiwvector3D,
            [](FlowField &flowField, int i, int j, int k, double &p,
               std::vector<int> & table) {
              *(&p + table[0]) = flowField.getPressure().getScalar(i, j) -
                                 flowField.getU(i, j) * flowField.getU(i, j);
              *(&p + table[1]) = flowField.getVelocity().getVector(i, j)[0];
              *(&p + table[2]) = flowField.getVelocity().getVector(i, j)[1];
              *(&p + table[3]) = flowField.getNu(i, j);
              *(&p + table[4]) = flowField.getPressure().getScalar(i, j);
            },
            [](FlowField &flowField, int i, int j, int k, double &p,
               std::vector<int> & table) {
              *(&p + table[0]) =
                  flowField.getPressure().getScalar(i, j, k) -
                  flowField.getU(i, j, k) * flowField.getU(i, j, k);
              *(&p + table[1]) = flowField.getVelocity().getVector(i, j, k)[0];
              *(&p + table[2]) = flowField.getVelocity().getVector(i, j, k)[1];
              *(&p + table[3]) = flowField.getVelocity().getVector(i, j, k)[2];
              *(&p + table[4]) = flowField.getNu(i, j, k);
              *(&p + table[5]) = flowField.getPressure().getScalar(i, j, k);
            }),
        _mpiir(*_flowField, parameters, _mpiiwvector2D, _mpiiwvector3D,
               [](FlowField &flowField, int i, int j, int k, double &p,
                  std::vector<int> & table) {
                 flowField.getPressure().getScalar(i, j) =
                     table[4] != -1 ? *(&p + table[4]) : 0.0;
                 flowField.getVelocity().getVector(i, j)[0] =
                     table[1] != -1 ? *(&p + table[1]) : 0.0;
                 flowField.getVelocity().getVector(i, j)[1] =
                     table[2] != -1 ? *(&p + table[2]) : 0.0;
                 flowField.getNu(i, j) =
                     table[3] != -1 ? *(&p + table[3]) : 0.0;
               },
               [](FlowField &flowField, int i, int j, int k, double &p,
                  std::vector<int> & table) {
                 flowField.getPressure().getScalar(i, j, k) =
                     table[5] != -1 ? *(&p + table[5]) : 0.0;
                 flowField.getVelocity().getVector(i, j, k)[0] =
                     table[1] != -1 ? *(&p + table[1]) : 0.0;
                 flowField.getVelocity().getVector(i, j, k)[1] =
                     table[2] != -1 ? *(&p + table[2]) : 0.0;
                 flowField.getVelocity().getVector(i, j, k)[2] =
                     table[3] != -1 ? *(&p + table[3]) : 0.0;
                 flowField.getNu(i, j, k) =
                     table[4] != -1 ? *(&p + table[4]) : 0.0;
               }),
        _nutst(parameters),
        _nutit(*_flowField, _parameters, _nutst, 1, 0),
        _hst(parameters, gm),
        _hit(*_flowField, _parameters, _hst, 0, 0),
        _maxmnutst(parameters),
        _maxmnutit(*_flowField, _parameters, _maxmnutst, 1, 0) {
    // distance to the next wall
    this->scalarStencils.push_back(CellDataStencil<double, FlowField>(
        this->_parameters, "h",
        [](FlowField &f, int i, int j) { return f.getH(i, j); },
        [](FlowField &f, int i, int j, int k) { return f.getH(i, j, k); }));

    // vortex viscosity
    this->scalarStencils.push_back(CellDataStencil<double, FlowField>(
        this->_parameters, "nu",
        [](FlowField &f, int i, int j) { return f.getNu(i, j); },
        [](FlowField &f, int i, int j, int k) { return f.getNu(i, j, k); }));

    // mixing length
    this->scalarStencils.push_back(CellDataStencil<double, FlowField>(
        this->_parameters, "lm",
        [](FlowField &f, int i, int j) { return f.getLm(i, j); },
        [](FlowField &f, int i, int j, int k) { return f.getLm(i, j, k); }));

    // velocity fluctuation
    this->scalarStencils.push_back(CellDataStencil<double, FlowField>(
        this->_parameters, "u",
        [](FlowField &f, int i, int j) { return f.getU(i, j); },
        [](FlowField &f, int i, int j, int k) { return f.getU(i, j, k); }));

    // actual pressure (without tke)
    this->scalarStencils.push_back(CellDataStencil<double, FlowField>(
        this->_parameters, "pressureactual",
        [](FlowField &f, int i, int j) {
          return f.getPressure().getScalar(i, j) - f.getU(i, j) * f.getU(i, j);
        },
        [](FlowField &f, int i, int j, int k) {
          return f.getPressure().getScalar(i, j, k) -
                 f.getU(i, j, k) * f.getU(i, j, k);
        }));
  }

  virtual ~Simulation() {}

  /** initialises the flow field according to the scenario */
  virtual void initializeFlowField() {
    if (_parameters.simulation.scenario == "taylor-green") {
      // currently, a particular initialization is only requrid for the
      // taylor-green vortex
      InitTaylorGreenFlowFieldStencil stencil(_parameters);
      FieldIterator<nseof::FlowField> iterator(*_flowField, _parameters,
                                               stencil);
      iterator.iterate();
    } else if (_parameters.simulation.scenario == "channel") {
      BFStepInitStencil stencil(_parameters);
      FieldIterator<nseof::FlowField> iterator(*_flowField, _parameters,
                                               stencil, 0, 1);
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
      FieldIterator<nseof::FlowField> iterator(*_flowField, _parameters,
                                               stencil, 0, 1);
      iterator.iterate();
    }
    _solver.reInitMatrix();

    _hit.iterate();
  }

  virtual void init() {
    // from super
    FlowFieldSimulation::init();

    // compute vortex viscosity
    _nutit.iterate();

    // communicate vortex viscosity
    nutComm.communicate(*this->_flowField);
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

  virtual void serialize() { _mpiiw.iterate(); }

  virtual void deserialize() { _mpiir.iterate(); }

 protected:
  /** sets the time step*/
  virtual void setTimeStep() {
    _maxmnutst.reset();
    _maxmnutit.iterate();

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
                                               _maxmnutst.getMaximum()) /
                                              (2 * factor),
                                          1.0 / _maxUStencil.getMaxValues()[0]),
                                 1.0 / _maxUStencil.getMaxValues()[1]));

    globalMin = MY_FLOAT_MAX;
    MPI_Allreduce(&localMin, &globalMin, 1, MY_MPI_FLOAT, MPI_MIN,
                  PETSC_COMM_WORLD);

    _parameters.timestep.dt = globalMin;
    _parameters.timestep.dt *= _parameters.timestep.tau;
    _parameters.timestep.dt =
        std::min(_parameters.timestep.dt, _parameters.timestep.dtu);
  }
};
}
}
}

#endif  // _NSEOF_FLOWMODELS_ALGEBRAIC_SIMULATION_H_
