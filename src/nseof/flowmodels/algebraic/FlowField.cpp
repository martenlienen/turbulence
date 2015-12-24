#include "FlowField.h"

namespace nseof {

namespace flowmodels {

namespace algebraic {

FlowField::FlowField(const Parameters& parameters)
    : nseof::FlowField(parameters),
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
