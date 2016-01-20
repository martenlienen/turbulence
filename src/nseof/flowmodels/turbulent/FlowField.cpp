// Ensure that mpi.h is included before everything else (otherwise MAC-Cluster
// complains)
#include "../../plotting/LambdaReader.h"

#include <memory>

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
              : ScalarField(getCellsX(), getCellsY(), getCellsZ())) {
  readers.push_back(
      std::make_unique<nseof::plotting::LambdaReader<FLOAT, FlowField, 1>>(
          *this, "h", [](FlowField& ff, int i, int j, int k) {
            return std::array<FLOAT, 1>{ff.getH(i, j, k)};
          }));

  readers.push_back(
      std::make_unique<nseof::plotting::LambdaReader<FLOAT, FlowField, 1>>(
          *this, "nu", [](FlowField& ff, int i, int j, int k) {
            return std::array<FLOAT, 1>{ff.getNu(i, j, k)};
          }));

  // mixing length
  readers.push_back(
      std::make_unique<nseof::plotting::LambdaReader<FLOAT, FlowField, 1>>(
          *this, "lm", [](FlowField& ff, int i, int j, int k) {
            return std::array<FLOAT, 1>{ff.getLm(i, j, k)};
          }));

  // velocity fluctuation
  readers.push_back(
      std::make_unique<nseof::plotting::LambdaReader<FLOAT, FlowField, 1>>(
          *this, "u", [](FlowField& ff, int i, int j, int k) {
            return std::array<FLOAT, 1>{ff.getU(i, j, k)};
          }));

  readers.push_back(
      std::make_unique<nseof::plotting::LambdaReader<FLOAT, FlowField, 1>>(
          *this, "Pressure-Actual", [](FlowField& ff, int i, int j, int k) {
            FLOAT u = ff.getU(i, j, k);
            FLOAT p = ff.getPressure().getScalar(i, j, k) - u * u;

            return std::array<FLOAT, 1>{p};
          }));
}

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
