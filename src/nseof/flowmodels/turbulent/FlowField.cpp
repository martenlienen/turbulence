#include "FlowField.h"

namespace nseof {

namespace flowmodels {

namespace turbulent {

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
}
}
}
