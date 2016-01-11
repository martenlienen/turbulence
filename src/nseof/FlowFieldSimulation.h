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
#include "Definitions.h"

#include "stencils/CellDataStencil.h"

#include "LinearSolver.h"
#include "solvers/SORSolver.h"
#include "solvers/PetscSolver.h"

#include "parallelManagers/MPICommunicator.h"

#include "MultiTimer.h"

#define GHOST_OFFSET 2

namespace nseof {

template <typename FF>
class FlowFieldSimulation : public Simulation {
 protected:
  std::unique_ptr<FF> _flowField;

  MPICommunicator<FLOAT, FF> pressureComm{
      *this->_flowField, this->_parameters,
      [](FF &flowField, int i, int j, int k, FLOAT &p) {
        p = flowField.getPressure().getScalar(i, j, k);
      },
      [](FF &flowField, int i, int j, int k, FLOAT &p) {
        flowField.getPressure().getScalar(i, j, k) = p;
      }};

  MPICommunicator<std::array<FLOAT, 3>, FF> velocityComm{
      *this->_flowField, this->_parameters,
      this->_parameters.geometry.dim==2?
        [](FF &flowField, int i, int j, int k, std::array<FLOAT, 3> &v) {
        std::copy_n(flowField.getVelocity().getVector(i, j, k),
                    2, v.data());
      }:[](FF &flowField, int i, int j, int k, std::array<FLOAT, 3> &v) {
        std::copy_n(flowField.getVelocity().getVector(i, j, k),
                    3, v.data());
      },
      this->_parameters.geometry.dim==2?[](FF &flowField, int i, int j, int k, std::array<FLOAT, 3> &v) {
        std::copy_n(v.data(), 2,
                    flowField.getVelocity().getVector(i, j, k));
      }:
        [](FF &flowField, int i, int j, int k, std::array<FLOAT, 3> &v) {
        std::copy_n(v.data(), 3,
                    flowField.getVelocity().getVector(i, j, k));
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
          }),
      CellDataStencil<std::vector<double>, FF>(
          this->_parameters, "velocity-raw",
          [](FF &f, int i, int j) {
            FLOAT *v = f.getVelocity().getVector(i, j);

            return std::vector<double>{v[0], v[1], 0};
          },
          [](FF &f, int i, int j, int k) {
            FLOAT *v = f.getVelocity().getVector(i, j, k);

            return std::vector<double>{v[0], v[1], v[2]};
          })};
          
          

 public:
  FlowFieldSimulation(Parameters &parameters, FF* flowField)
      : Simulation(parameters), _flowField(flowField) {}

  virtual ~FlowFieldSimulation() {}

  /** TODO WS1: plots the flow field. */
  virtual void plotVTK(int rank, int timeStep) {
    if (!this->_parameters.vtk.enabled) {
      return;
    }

    // TODO WS1: create VTKStencil and respective iterator; iterate stencil
    //           over _flowField and write flow field information to vtk file

    int los = _parameters.vtk.lowoffset;
    int hos = _parameters.vtk.highoffset;

    vtk::Dataset dataset = this->datasetFromMesh(los, hos);
    std::vector<std::unique_ptr<vtk::CellData>> cellData;
    cellData.reserve(this->scalarStencils.size() + this->vectorStencils.size());

    for (auto &stencil : this->scalarStencils) {
      stencil.reset();
      FieldIterator<FF> iterator(*this->_flowField, this->_parameters, stencil,
                                 los - 1, 1 - hos);
      iterator.iterate();
      cellData.push_back(stencil.get());
    }

    for (auto &stencil : this->vectorStencils) {
      stencil.reset();
      FieldIterator<FF> iterator(*this->_flowField, this->_parameters, stencil,
                                 los - 1, 1 - hos);
      iterator.iterate();
      cellData.push_back(stencil.get());
    }

    vtk::File file(dataset, std::move(cellData));
    file.write(this->_parameters.vtk.prefix, rank, timeStep);
  }
  
  virtual void serialize(){}
  virtual void deserialize(){}

 private:
  vtk::Dataset datasetFromMesh(int los, int hos) {
    Parameters &p = this->_parameters;
    ParallelParameters &pp = p.parallel;
    GeometricParameters &gp = p.geometry;
    Meshsize *mesh = p.meshsize;

    int cellsX = pp.localSize[0] - (los + hos - 3);
    int cellsY = pp.localSize[1] - (los + hos - 3);
    int cellsZ = gp.dim == 3 ? pp.localSize[2] - (los + hos - 3) : 1;
    int pointsX = cellsX + 1;
    int pointsY = cellsY + 1;
    int pointsZ = gp.dim == 3 ? cellsZ + 1 : 1;
    int numPoints = pointsX * pointsY * pointsZ;

    std::vector<vtk::Point> points(numPoints);

    int losz = _parameters.geometry.dim == 2 ? 0 : los;

    for (int k = losz, p = 0; k < pointsZ + losz; k++) {
      for (int j = los; j < pointsY + los; j++) {
        for (int i = los; i < pointsX + los; i++, p++) {
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
}

#endif  // _FLOW_FIELD_SIMULATION_H_
