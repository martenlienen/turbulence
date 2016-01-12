#include "NeumannBoundaryStencils.h"

namespace nseof {

NeumannVelocityBoundaryStencil::NeumannVelocityBoundaryStencil(
    const Parameters& parameters)
    : BoundaryStencil<FlowField>(parameters) {}

void NeumannVelocityBoundaryStencil::applyLeftWall(FlowField& flowField, int i,
                                                   int j) {
  flowField.getVelocity().getVector(i - 1, j)[0] =
      flowField.getVelocity().getVector(i, j)[0];
  flowField.getVelocity().getVector(i, j)[1] =
      flowField.getVelocity().getVector(i + 1, j)[1];
}

void NeumannVelocityBoundaryStencil::applyRightWall(FlowField& flowField, int i,
                                                    int j) {
  flowField.getVelocity().getVector(i, j)[0] =
      flowField.getVelocity().getVector(i - 1, j)[0];
  flowField.getVelocity().getVector(i, j)[1] =
      flowField.getVelocity().getVector(i - 1, j)[1];
}

void NeumannVelocityBoundaryStencil::applyBottomWall(FlowField& flowField,
                                                     int i, int j) {
  flowField.getVelocity().getVector(i, j)[0] =
      flowField.getVelocity().getVector(i, j + 1)[0];
  flowField.getVelocity().getVector(i, j - 1)[1] =
      flowField.getVelocity().getVector(i, j)[1];
}

void NeumannVelocityBoundaryStencil::applyTopWall(FlowField& flowField, int i,
                                                  int j) {
  flowField.getVelocity().getVector(i, j)[0] =
      flowField.getVelocity().getVector(i, j - 1)[0];
  flowField.getVelocity().getVector(i, j)[1] =
      flowField.getVelocity().getVector(i, j - 1)[1];
}

void NeumannVelocityBoundaryStencil::applyLeftWall(FlowField& flowField, int i,
                                                   int j, int k) {
  flowField.getVelocity().getVector(i - 1, j, k)[0] =
      flowField.getVelocity().getVector(i, j, k)[0];
  flowField.getVelocity().getVector(i, j, k)[1] =
      flowField.getVelocity().getVector(i + 1, j, k)[1];
  flowField.getVelocity().getVector(i, j, k)[2] =
      flowField.getVelocity().getVector(i + 1, j, k)[2];
}

void NeumannVelocityBoundaryStencil::applyRightWall(FlowField& flowField, int i,
                                                    int j, int k) {
  flowField.getVelocity().getVector(i, j, k)[0] =
      flowField.getVelocity().getVector(i - 1, j, k)[0];
  flowField.getVelocity().getVector(i, j, k)[1] =
      flowField.getVelocity().getVector(i - 1, j, k)[1];
  flowField.getVelocity().getVector(i, j, k)[2] =
      flowField.getVelocity().getVector(i - 1, j, k)[2];
}

void NeumannVelocityBoundaryStencil::applyBottomWall(FlowField& flowField,
                                                     int i, int j, int k) {
  flowField.getVelocity().getVector(i, j, k)[0] =
      flowField.getVelocity().getVector(i, j + 1, k)[0];
  flowField.getVelocity().getVector(i, j - 1, k)[1] =
      flowField.getVelocity().getVector(i, j, k)[1];
  flowField.getVelocity().getVector(i, j, k)[2] =
      flowField.getVelocity().getVector(i, j + 1, k)[2];
}

void NeumannVelocityBoundaryStencil::applyTopWall(FlowField& flowField, int i,
                                                  int j, int k) {
  flowField.getVelocity().getVector(i, j, k)[0] =
      flowField.getVelocity().getVector(i, j - 1, k)[0];
  flowField.getVelocity().getVector(i, j, k)[1] =
      flowField.getVelocity().getVector(i, j - 1, k)[1];
  flowField.getVelocity().getVector(i, j, k)[2] =
      flowField.getVelocity().getVector(i, j - 1, k)[2];
}

void NeumannVelocityBoundaryStencil::applyFrontWall(FlowField& flowField, int i,
                                                    int j, int k) {
  flowField.getVelocity().getVector(i, j, k)[0] =
      flowField.getVelocity().getVector(i, j, k + 1)[0];
  flowField.getVelocity().getVector(i, j, k)[1] =
      flowField.getVelocity().getVector(i, j, k + 1)[1];
  flowField.getVelocity().getVector(i, j, k - 1)[2] =
      flowField.getVelocity().getVector(i, j, k)[2];
}

void NeumannVelocityBoundaryStencil::applyBackWall(FlowField& flowField, int i,
                                                   int j, int k) {
  flowField.getVelocity().getVector(i, j, k)[0] =
      flowField.getVelocity().getVector(i, j, k - 1)[0];
  flowField.getVelocity().getVector(i, j, k)[1] =
      flowField.getVelocity().getVector(i, j, k - 1)[1];
  flowField.getVelocity().getVector(i, j, k)[2] =
      flowField.getVelocity().getVector(i, j, k - 1)[2];
}

NeumannFGHBoundaryStencil::NeumannFGHBoundaryStencil(
    const Parameters& parameters)
    : BoundaryStencil<FlowField>(parameters) {}

// These are left empty. The right values should be computed by the FGH body
// stencil

void NeumannFGHBoundaryStencil::applyLeftWall(FlowField& flowField, int i,
                                              int j) {}
void NeumannFGHBoundaryStencil::applyRightWall(FlowField& flowField, int i,
                                               int j) {}
void NeumannFGHBoundaryStencil::applyBottomWall(FlowField& flowField, int i,
                                                int j) {}
void NeumannFGHBoundaryStencil::applyTopWall(FlowField& flowField, int i,
                                             int j) {}

// 3D stencils

void NeumannFGHBoundaryStencil::applyLeftWall(FlowField& flowField, int i,
                                              int j, int k) {}
void NeumannFGHBoundaryStencil::applyRightWall(FlowField& flowField, int i,
                                               int j, int k) {}
void NeumannFGHBoundaryStencil::applyBottomWall(FlowField& flowField, int i,
                                                int j, int k) {}
void NeumannFGHBoundaryStencil::applyTopWall(FlowField& flowField, int i, int j,
                                             int k) {}
void NeumannFGHBoundaryStencil::applyFrontWall(FlowField& flowField, int i,
                                               int j, int k) {}
void NeumannFGHBoundaryStencil::applyBackWall(FlowField& flowField, int i,
                                              int j, int k) {}

SymmetryVelocityBoundaryStencil::SymmetryVelocityBoundaryStencil(
    const Parameters& parameters)
    : BoundaryStencil<FlowField>(parameters) {}

void SymmetryVelocityBoundaryStencil::applyLeftWall(FlowField& flowField, int i,
                                                    int j) {}

void SymmetryVelocityBoundaryStencil::applyRightWall(FlowField& flowField,
                                                     int i, int j) {}

void SymmetryVelocityBoundaryStencil::applyBottomWall(FlowField& flowField,
                                                      int i, int j) {}

void SymmetryVelocityBoundaryStencil::applyTopWall(FlowField& flowField, int i,
                                                   int j) {
  flowField.getVelocity().getVector(i, j)[0] =
      flowField.getVelocity().getVector(i, j - 1)[0];
  flowField.getVelocity().getVector(i, j)[1] =
      -flowField.getVelocity().getVector(i, j - 2)[1];
}

void SymmetryVelocityBoundaryStencil::applyLeftWall(FlowField& flowField, int i,
                                                    int j, int k) {}

void SymmetryVelocityBoundaryStencil::applyRightWall(FlowField& flowField,
                                                     int i, int j, int k) {}

void SymmetryVelocityBoundaryStencil::applyBottomWall(FlowField& flowField,
                                                      int i, int j, int k) {}

void SymmetryVelocityBoundaryStencil::applyTopWall(FlowField& flowField, int i,
                                                   int j, int k) {
  flowField.getVelocity().getVector(i, j, k)[0] =
      flowField.getVelocity().getVector(i, j - 1, k)[0];
  flowField.getVelocity().getVector(i, j, k)[1] =
      -flowField.getVelocity().getVector(i, j - 2, k)[1];
  flowField.getVelocity().getVector(i, j, k)[2] =
      +flowField.getVelocity().getVector(i, j - 2, k)[2];
}

void SymmetryVelocityBoundaryStencil::applyFrontWall(FlowField& flowField,
                                                     int i, int j, int k) {}

void SymmetryVelocityBoundaryStencil::applyBackWall(FlowField& flowField, int i,
                                                    int j, int k) {
  flowField.getVelocity().getVector(i, j, k)[0] =
      flowField.getVelocity().getVector(i, j, k - 1)[0];
  flowField.getVelocity().getVector(i, j, k)[1] =
      +flowField.getVelocity().getVector(i, j, k - 2)[1];
  flowField.getVelocity().getVector(i, j, k)[2] =
      -flowField.getVelocity().getVector(i, j, k - 2)[2];
}

SymmetryFGHBoundaryStencil::SymmetryFGHBoundaryStencil(
    const Parameters& parameters)
    : BoundaryStencil<FlowField>(parameters) {}

// These are left empty. The right values should be computed by the FGH body
// stencil

void SymmetryFGHBoundaryStencil::applyLeftWall(FlowField& flowField, int i,
                                               int j) {}
void SymmetryFGHBoundaryStencil::applyRightWall(FlowField& flowField, int i,
                                                int j) {}
void SymmetryFGHBoundaryStencil::applyBottomWall(FlowField& flowField, int i,
                                                 int j) {}
void SymmetryFGHBoundaryStencil::applyTopWall(FlowField& flowField, int i,
                                              int j) {
  flowField.getVelocity().getVector(i, j)[0] =
      flowField.getVelocity().getVector(i, j - 1)[0];
  flowField.getVelocity().getVector(i, j - 1)[1] = 0;
}

void SymmetryFGHBoundaryStencil::applyLeftWall(FlowField& flowField, int i,
                                               int j, int k) {}
void SymmetryFGHBoundaryStencil::applyRightWall(FlowField& flowField, int i,
                                                int j, int k) {}
void SymmetryFGHBoundaryStencil::applyBottomWall(FlowField& flowField, int i,
                                                 int j, int k) {}
void SymmetryFGHBoundaryStencil::applyTopWall(FlowField& flowField, int i,
                                              int j, int k) {
  flowField.getVelocity().getVector(i, j, k)[0] =
      flowField.getVelocity().getVector(i, j - 1, k)[0];
  flowField.getVelocity().getVector(i, j - 1, k)[1] = 0;
  flowField.getVelocity().getVector(i, j, k)[0] =
      flowField.getVelocity().getVector(i, j - 1, k)[2];
}
void SymmetryFGHBoundaryStencil::applyFrontWall(FlowField& flowField, int i,
                                                int j, int k) {}
void SymmetryFGHBoundaryStencil::applyBackWall(FlowField& flowField, int i,
                                               int j, int k) {
  flowField.getVelocity().getVector(i, j, k)[0] =
      flowField.getVelocity().getVector(i, j, k - 1)[0];
  flowField.getVelocity().getVector(i, j, k - 1)[1] = 0;
  flowField.getVelocity().getVector(i, j, k)[0] =
      flowField.getVelocity().getVector(i, j, k - 1)[2];
}
}
