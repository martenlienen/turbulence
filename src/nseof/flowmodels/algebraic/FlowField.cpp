#include "FlowField.h"

namespace nseof {

namespace flowmodels {

namespace algebraic {

FlowField::FlowField(const Parameters& parameters)
    : _flowField(parameters),
      _nu(parameters.geometry.dim == 2
              ? ScalarField(getCellsX(), getCellsY())
              : ScalarField(getCellsX(), getCellsY(), getCellsZ())),
      _h(parameters.geometry.dim == 2
             ? ScalarField(getCellsX(), getCellsY())
             : ScalarField(getCellsX(), getCellsY(), getCellsZ())),
      _u(parameters.geometry.dim == 2
             ? ScalarField(getCellsX(), getCellsY())
             : ScalarField(getCellsX(), getCellsY(), getCellsZ())),
      _lm(parameters.geometry.dim == 2
              ? ScalarField(getCellsX(), getCellsY())
              : ScalarField(getCellsX(), getCellsY(), getCellsZ())) {}

int FlowField::getNx() const { return _flowField.getNx(); }

int FlowField::getNy() const { return _flowField.getNy(); }

int FlowField::getNz() const { return _flowField.getNz(); }

int FlowField::getCellsX() const { return _flowField.getCellsX(); }

int FlowField::getCellsY() const { return _flowField.getCellsY(); }

int FlowField::getCellsZ() const { return _flowField.getCellsZ(); }

ScalarField& FlowField::getPressure() { return _flowField.getPressure(); }

VectorField& FlowField::getVelocity() { return _flowField.getVelocity(); }

IntScalarField& FlowField::getFlags() { return _flowField.getFlags(); }

VectorField& FlowField::getFGH() { return _flowField.getFGH(); }

ScalarField& FlowField::getRHS() { return _flowField.getRHS(); }

void FlowField::getPressureAndVelocity(FLOAT& pressure, FLOAT* const velocity,
                                       int i, int j) {
  _flowField.getPressureAndVelocity(pressure, velocity, i, j);
}

void FlowField::getPressureAndVelocity(FLOAT& pressure, FLOAT* const velocity,
                                       int i, int j, int k) {
  _flowField.getPressureAndVelocity(pressure, velocity, i, j, k);
}

FLOAT& FlowField::getNu(int i, int j) { return _nu.getScalar(i, j); }

FLOAT& FlowField::getNu(int i, int j, int k) { return _nu.getScalar(i, j, k); }

FLOAT& FlowField::getH(int i, int j) { return _h.getScalar(i, j); }

FLOAT& FlowField::getH(int i, int j, int k) { return _h.getScalar(i, j, k); }

FLOAT& FlowField::getLm(int i, int j) { return _lm.getScalar(i, j); }

FLOAT& FlowField::getLm(int i, int j, int k) { return _lm.getScalar(i, j, k); }

FLOAT& FlowField::getU(int i, int j) { return _u.getScalar(i, j); }
FLOAT& FlowField::getU(int i, int j, int k) { return _u.getScalar(i, j, k); }

void loadLocalNu2D(const Parameters& parameters, FlowField& flowField,
                   FLOAT* const localNu, int i, int j) {
  const FLOAT nu = 1 / parameters.flow.Re;
  for (int row = -1; row <= 1; row++) {
    for (int column = -1; column <= 1; column++) {
      FLOAT nut = flowField.getNu(i + column, j + row);
      localNu[39 + 9 * row + 3 * column] = nut + nu;
    }
  }
}

void loadLocalNu3D(const Parameters& parameters, FlowField& flowField,
                   FLOAT* const localNu, int i, int j, int k) {
  const FLOAT nu = 1 / parameters.flow.Re;

  for (int layer = -1; layer <= 1; layer++) {
    for (int row = -1; row <= 1; row++) {
      for (int column = -1; column <= 1; column++) {
        FLOAT nut = flowField.getNu(i + column, j + row, k + layer);
        localNu[39 + 27 * layer + 9 * row + 3 * column] = nut + nu;
      }
    }
  }
}
}
}
}
