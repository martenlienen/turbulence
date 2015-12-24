#include "FlowField.h"
#include "FlowFieldTurbA.h"

namespace nseof {

FlowFieldTurbA::FlowFieldTurbA(const Parameters& parameters)
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

int FlowFieldTurbA::getNx() const { return _flowField.getNx(); }

int FlowFieldTurbA::getNy() const { return _flowField.getNy(); }

int FlowFieldTurbA::getNz() const { return _flowField.getNz(); }

int FlowFieldTurbA::getCellsX() const { return _flowField.getCellsX(); }

int FlowFieldTurbA::getCellsY() const { return _flowField.getCellsY(); }

int FlowFieldTurbA::getCellsZ() const { return _flowField.getCellsZ(); }

ScalarField& FlowFieldTurbA::getPressure() { return _flowField.getPressure(); }

VectorField& FlowFieldTurbA::getVelocity() { return _flowField.getVelocity(); }

IntScalarField& FlowFieldTurbA::getFlags() { return _flowField.getFlags(); }

VectorField& FlowFieldTurbA::getFGH() { return _flowField.getFGH(); }

ScalarField& FlowFieldTurbA::getRHS() { return _flowField.getRHS(); }

void FlowFieldTurbA::getPressureAndVelocity(FLOAT& pressure,
                                            FLOAT* const velocity, int i,
                                            int j) {
  _flowField.getPressureAndVelocity(pressure, velocity, i, j);
}

void FlowFieldTurbA::getPressureAndVelocity(FLOAT& pressure,
                                            FLOAT* const velocity, int i, int j,
                                            int k) {
  _flowField.getPressureAndVelocity(pressure, velocity, i, j, k);
}

FLOAT& FlowFieldTurbA::getNu(int i, int j) { return _nu.getScalar(i, j); }

FLOAT& FlowFieldTurbA::getNu(int i, int j, int k) {
  return _nu.getScalar(i, j, k);
}

FLOAT& FlowFieldTurbA::getH(int i, int j) { return _h.getScalar(i, j); }

FLOAT& FlowFieldTurbA::getH(int i, int j, int k) {
  return _h.getScalar(i, j, k);
}

FLOAT& FlowFieldTurbA::getLm(int i, int j) { return _lm.getScalar(i, j); }

FLOAT& FlowFieldTurbA::getLm(int i, int j, int k) {
  return _lm.getScalar(i, j, k);
}

FLOAT& FlowFieldTurbA::getU(int i, int j) { return _u.getScalar(i, j); }
FLOAT& FlowFieldTurbA::getU(int i, int j, int k) {
  return _u.getScalar(i, j, k);
}
}
