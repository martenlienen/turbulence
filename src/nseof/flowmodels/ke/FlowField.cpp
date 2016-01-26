// Ensure that mpi.h is included before everything else (otherwise MAC-Cluster
// complains)
#include "../../plotting/LambdaReader.h"

#include <memory>

#include "FlowField.h"

namespace nseof {

namespace flowmodels {

namespace ke {

FlowField::FlowField(const Parameters& parameters)
    : nseof::flowmodels::turbulent::FlowField(parameters),
      _fmu(parameters.geometry.dim == 2
               ? ScalarField(this->getCellsX(), this->getCellsY())
               : ScalarField(this->getCellsX(), this->getCellsY(),
                             this->getCellsZ())),
      _f1(parameters.geometry.dim == 2
              ? ScalarField(this->getCellsX(), this->getCellsY())
              : ScalarField(this->getCellsX(), this->getCellsY(),
                            this->getCellsZ())),
      _f2(parameters.geometry.dim == 2
              ? ScalarField(this->getCellsX(), this->getCellsY())
              : ScalarField(this->getCellsX(), this->getCellsY(),
                            this->getCellsZ())),
      _f3(parameters.geometry.dim == 2
              ? ScalarField(this->getCellsX(), this->getCellsY())
              : ScalarField(this->getCellsX(), this->getCellsY(),
                            this->getCellsZ())),
      _D(parameters.geometry.dim == 2
             ? ScalarField(this->getCellsX(), this->getCellsY())
             : ScalarField(this->getCellsX(), this->getCellsY(),
                           this->getCellsZ())),
      _E(parameters.geometry.dim == 2
             ? ScalarField(this->getCellsX(), this->getCellsY())
             : ScalarField(this->getCellsX(), this->getCellsY(),
                           this->getCellsZ())),
      _epsilon(parameters.geometry.dim == 2
                   ? ScalarField(this->getCellsX(), this->getCellsY())
                   : ScalarField(this->getCellsX(), this->getCellsY(),
                                 this->getCellsZ())),
      _tke(parameters.geometry.dim == 2
               ? ScalarField(this->getCellsX(), this->getCellsY())
               : ScalarField(this->getCellsX(), this->getCellsY(),
                             this->getCellsZ())),
      _rhs_epsilon(parameters.geometry.dim == 2
                       ? ScalarField(this->getCellsX(), this->getCellsY())
                       : ScalarField(this->getCellsX(), this->getCellsY(),
                                     this->getCellsZ())),
      _rhs_tke(parameters.geometry.dim == 2
                   ? ScalarField(this->getCellsX(), this->getCellsY())
                   : ScalarField(this->getCellsX(), this->getCellsY(),
                                 this->getCellsZ())),
      _sijsij(parameters.geometry.dim == 2
                  ? ScalarField(this->getCellsX(), this->getCellsY())
                  : ScalarField(this->getCellsX(), this->getCellsY(),
                                this->getCellsZ())) {
  // turbulent kinetic energy
  readers.push_back(
      std::make_unique<nseof::plotting::LambdaReader<FLOAT, FlowField, 1>>(
          *this, "tke", [](FlowField& ff, int i, int j, int k) {
            return std::array<FLOAT, 1>{ff.getTke(i, j, k)};
          }));

  // dissipation rate
  readers.push_back(
      std::make_unique<nseof::plotting::LambdaReader<FLOAT, FlowField, 1>>(
          *this, "epsilon", [](FlowField& ff, int i, int j, int k) {
            return std::array<FLOAT, 1>{ff.getEpsilon(i, j, k)};
          }));

  readers.push_back(
      std::make_unique<nseof::plotting::LambdaReader<FLOAT, FlowField, 1>>(
          *this, "F1", [](FlowField& ff, int i, int j, int k) {
            return std::array<FLOAT, 1>{ff.getF1(i, j, k)};
          }));

  readers.push_back(
      std::make_unique<nseof::plotting::LambdaReader<FLOAT, FlowField, 1>>(
          *this, "F2", [](FlowField& ff, int i, int j, int k) {
            return std::array<FLOAT, 1>{ff.getF2(i, j, k)};
          }));

  readers.push_back(
      std::make_unique<nseof::plotting::LambdaReader<FLOAT, FlowField, 1>>(
          *this, "Fmu", [](FlowField& ff, int i, int j, int k) {
            return std::array<FLOAT, 1>{ff.getFmu(i, j, k)};
          }));

  readers.push_back(
      std::make_unique<nseof::plotting::LambdaReader<FLOAT, FlowField, 1>>(
          *this, "SijSij", [](FlowField& ff, int i, int j, int k) {
            return std::array<FLOAT, 1>{ff.getsijsij(i, j, k)};
          }));

  readers.push_back(
      std::make_unique<nseof::plotting::LambdaReader<FLOAT, FlowField, 1>>(
          *this, "D", [](FlowField& ff, int i, int j, int k) {
            return std::array<FLOAT, 1>{ff.getD(i, j, k)};
          }));

  readers.push_back(
      std::make_unique<nseof::plotting::LambdaReader<FLOAT, FlowField, 1>>(
          *this, "E", [](FlowField& ff, int i, int j, int k) {
            return std::array<FLOAT, 1>{ff.getE(i, j, k)};
          }));
}

ScalarField& FlowField::getFmu() { return _fmu; }
ScalarField& FlowField::getF1() { return _f1; }
ScalarField& FlowField::getF2() { return _f2; }
ScalarField& FlowField::getF3() { return _f3; }
ScalarField& FlowField::getD() { return _D; }
ScalarField& FlowField::getE() { return _E; }
ScalarField& FlowField::getEpsilon() { return _epsilon; }
ScalarField& FlowField::getTke() { return _tke; }
ScalarField& FlowField::getRHSEpsilon() { return _rhs_epsilon; }
ScalarField& FlowField::getRHSTke() { return _rhs_tke; }
ScalarField& FlowField::getsijsij() { return _sijsij; }

FLOAT& FlowField::getFmu(int i, int j) { return _fmu.getScalar(i, j); }

FLOAT& FlowField::getFmu(int i, int j, int k) {
  return _fmu.getScalar(i, j, k);
}
FLOAT& FlowField::getF1(int i, int j) { return _f1.getScalar(i, j); }

FLOAT& FlowField::getF1(int i, int j, int k) { return _f1.getScalar(i, j, k); }
FLOAT& FlowField::getF2(int i, int j) { return _f2.getScalar(i, j); }

FLOAT& FlowField::getF2(int i, int j, int k) { return _f2.getScalar(i, j, k); }
FLOAT& FlowField::getF3(int i, int j) { return _f3.getScalar(i, j); }

FLOAT& FlowField::getF3(int i, int j, int k) { return _f3.getScalar(i, j, k); }
FLOAT& FlowField::getD(int i, int j) { return _D.getScalar(i, j); }

FLOAT& FlowField::getD(int i, int j, int k) { return _D.getScalar(i, j, k); }
FLOAT& FlowField::getE(int i, int j) { return _E.getScalar(i, j); }

FLOAT& FlowField::getE(int i, int j, int k) { return _E.getScalar(i, j, k); }
FLOAT& FlowField::getEpsilon(int i, int j) { return _epsilon.getScalar(i, j); }

FLOAT& FlowField::getEpsilon(int i, int j, int k) {
  return _epsilon.getScalar(i, j, k);
}
FLOAT& FlowField::getTke(int i, int j) { return _tke.getScalar(i, j); }

FLOAT& FlowField::getTke(int i, int j, int k) {
  return _tke.getScalar(i, j, k);
}

FLOAT& FlowField::getRHSEpsilon(int i, int j) {
  return _rhs_epsilon.getScalar(i, j);
}

FLOAT& FlowField::getRHSEpsilon(int i, int j, int k) {
  return _rhs_epsilon.getScalar(i, j, k);
}
FLOAT& FlowField::getRHSTke(int i, int j) { return _rhs_tke.getScalar(i, j); }

FLOAT& FlowField::getRHSTke(int i, int j, int k) {
  return _rhs_tke.getScalar(i, j, k);
}

FLOAT& FlowField::getsijsij(int i, int j) { return _sijsij.getScalar(i, j); }

FLOAT& FlowField::getsijsij(int i, int j, int k) {
  return _rhs_tke.getScalar(i, j, k);
}
}
}
}
