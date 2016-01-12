#include "BFInputStencils.h"

namespace nseof {

namespace flowmodels {

namespace ke {

BFInputStencil::BFInputStencil(const Parameters& parameters)
    : BoundaryStencil<FlowField>(parameters) {}

// Most of the functions are empty, and they shouldn't be called, assuming that
// the input is always
// located at the left.

void BFInputStencil::applyLeftWall(FlowField& flowField, int i, int j) {
  FLOAT Uin = _parameters.walls.vectorLeft[0];
  FLOAT kin = 0.003 * Uin * Uin;
  FLOAT ein = _parameters.kEpsilon.cmu * pow(kin, 1.5) /
              (0.03 * _parameters.geometry.lengthY);

  flowField.getTke(i, j) = 2 * kin - flowField.getTke(i + 1, j);
  flowField.getEpsilon(i, j) = 2 * ein - flowField.getEpsilon(i + 1, j);
}
void BFInputStencil::applyRightWall(FlowField& flowField, int i, int j) {
  flowField.getTke(i, j) = flowField.getTke(i - 1, j);
  flowField.getEpsilon(i, j) = flowField.getEpsilon(i - 1, j);
}
void BFInputStencil::applyBottomWall(FlowField& flowField, int i, int j) {
  flowField.getTke(i, j) = -flowField.getTke(i, j + 1);
  flowField.getEpsilon(i, j) = flowField.getEpsilon(i, j + 1);
}
void BFInputStencil::applyTopWall(FlowField& flowField, int i, int j) {
  flowField.getTke(i, j) = -flowField.getTke(i, j - 1);
  flowField.getEpsilon(i, j) = flowField.getEpsilon(i, j - 1);
}

void BFInputStencil::applyLeftWall(FlowField& flowField, int i, int j, int k) {
  FLOAT Uin = _parameters.walls.vectorLeft[0];
  FLOAT kin = 0.003 * Uin * Uin;
  FLOAT ein = _parameters.kEpsilon.cmu * pow(kin, 1.5) /
              (0.03 * _parameters.geometry.lengthY);

  flowField.getTke(i, j, k) = 2 * kin - flowField.getTke(i + 1, j, k);
  flowField.getEpsilon(i, j, k) = 2 * ein - flowField.getEpsilon(i + 1, j, k);
}
void BFInputStencil::applyRightWall(FlowField& flowField, int i, int j, int k) {
  flowField.getTke(i, j, k) = flowField.getTke(i - 1, j, k);
  flowField.getEpsilon(i, j, k) = flowField.getEpsilon(i - 1, j, k);
}
void BFInputStencil::applyBottomWall(FlowField& flowField, int i, int j,
                                     int k) {
  flowField.getTke(i, j, k) = -flowField.getTke(i, j + 1, k);
  flowField.getEpsilon(i, j, k) = flowField.getEpsilon(i, j + 1, k);
}
void BFInputStencil::applyTopWall(FlowField& flowField, int i, int j, int k) {
  flowField.getTke(i, j, k) = -flowField.getTke(i, j - 1, k);
  flowField.getEpsilon(i, j, k) = flowField.getEpsilon(i, j - 1, k);
}
void BFInputStencil::applyFrontWall(FlowField& flowField, int i, int j, int k) {
  flowField.getTke(i, j, k) = -flowField.getTke(i, j, k + 1);
  flowField.getEpsilon(i, j, k) = flowField.getEpsilon(i, j, k + 1);
}
void BFInputStencil::applyBackWall(FlowField& flowField, int i, int j, int k) {
  flowField.getTke(i, j, k) = -flowField.getTke(i, j - 1, k);
  flowField.getEpsilon(i, j, k) = flowField.getEpsilon(i, j - 1, k);
}

SymmStencil::SymmStencil(const Parameters& parameters)
    : BoundaryStencil<FlowField>(parameters) {}

void SymmStencil::applyLeftWall(FlowField& flowField, int i, int j) {}

void SymmStencil::applyRightWall(FlowField& flowField, int i, int j) {}

void SymmStencil::applyBottomWall(FlowField& flowField, int i, int j) {}

void SymmStencil::applyTopWall(FlowField& flowField, int i, int j) {
  flowField.getTke(i, j) = flowField.getTke(i, j - 1);
  flowField.getEpsilon(i, j) = flowField.getEpsilon(i, j - 1);
}

void SymmStencil::applyLeftWall(FlowField& flowField, int i, int j, int k) {}

void SymmStencil::applyRightWall(FlowField& flowField, int i, int j, int k) {}

void SymmStencil::applyBottomWall(FlowField& flowField, int i, int j, int k) {}

void SymmStencil::applyTopWall(FlowField& flowField, int i, int j, int k) {
  // TODO
}

void SymmStencil::applyFrontWall(FlowField& flowField, int i, int j, int k) {}

void SymmStencil::applyBackWall(FlowField& flowField, int i, int j, int k) {
  // TODO
}

WallStencil::WallStencil(const Parameters& parameters)
    : BoundaryStencil<FlowField>(parameters) {}

// Most of the functions are empty, and they shouldn't be called, assuming that
// the input is always
// located at the left.

void WallStencil::applyLeftWall(FlowField& flowField, int i, int j) {
  flowField.getTke(i, j) = -flowField.getTke(i + 1, j);
  flowField.getEpsilon(i, j) = flowField.getEpsilon(i + 1, j);
}
void WallStencil::applyRightWall(FlowField& flowField, int i, int j) {
  flowField.getTke(i, j) = -flowField.getTke(i - 1, j);
  flowField.getEpsilon(i, j) = flowField.getEpsilon(i - 1, j);
}
void WallStencil::applyBottomWall(FlowField& flowField, int i, int j) {
  flowField.getTke(i, j) = -flowField.getTke(i, j + 1);
  flowField.getEpsilon(i, j) = flowField.getEpsilon(i, j + 1);
}
void WallStencil::applyTopWall(FlowField& flowField, int i, int j) {
  flowField.getTke(i, j) = -flowField.getTke(i, j - 1);
  flowField.getEpsilon(i, j) = flowField.getEpsilon(i, j - 1);
}

void WallStencil::applyLeftWall(FlowField& flowField, int i, int j, int k) {
  flowField.getTke(i, j, k) = -flowField.getTke(i + 1, j, k);
  flowField.getEpsilon(i, j, k) = flowField.getEpsilon(i + 1, j, k);
}
void WallStencil::applyRightWall(FlowField& flowField, int i, int j, int k) {
  flowField.getTke(i, j, k) = -flowField.getTke(i - 1, j, k);
  flowField.getEpsilon(i, j, k) = flowField.getEpsilon(i - 1, j, k);
}
void WallStencil::applyBottomWall(FlowField& flowField, int i, int j, int k) {
  flowField.getTke(i, j, k) = -flowField.getTke(i, j + 1, k);
  flowField.getEpsilon(i, j, k) = flowField.getEpsilon(i, j + 1, k);
}
void WallStencil::applyTopWall(FlowField& flowField, int i, int j, int k) {
  flowField.getTke(i, j, k) = -flowField.getTke(i, j - 1, k);
  flowField.getEpsilon(i, j, k) = flowField.getEpsilon(i, j - 1, k);
}
void WallStencil::applyFrontWall(FlowField& flowField, int i, int j, int k) {
  flowField.getTke(i, j, k) = -flowField.getTke(i, j, k + 1);
  flowField.getEpsilon(i, j, k) = flowField.getEpsilon(i, j, k + 1);
}
void WallStencil::applyBackWall(FlowField& flowField, int i, int j, int k) {
  flowField.getTke(i, j, k) = -flowField.getTke(i, j - 1, k);
  flowField.getEpsilon(i, j, k) = flowField.getEpsilon(i, j - 1, k);
}
}
}
}
