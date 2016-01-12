#include <algorithm>

#include "../../Definitions.h"
#include "../../stencils/StencilFunctions.h"

#include "NutStencil.h"

namespace nseof {

namespace flowmodels {

namespace algebraic {

NutStencil::NutStencil(const Parameters& parameters)
    : FieldStencil<FlowField>(parameters) {
  // do you want to calculate vortex viscosity?
  if (parameters.flow.type == "laminar") {
    _fs = new NutStencilL(parameters);
  } else {
    _fs = new NutStencilA(parameters);
  }
}

NutStencil::~NutStencil() { delete _fs; }

void NutStencil::apply(FlowField& flowField, int i, int j) {
  _fs->apply(flowField, i, j);
}

void NutStencil::apply(FlowField& flowField, int i, int j, int k) {
  _fs->apply(flowField, i, j, k);
}

NutStencilL::NutStencilL(const Parameters& parameters)
    : FieldStencil<FlowField>(parameters) {}

void NutStencilL::apply(FlowField& flowField, int i, int j) {
  // no vortex viscosity, velociy fluctuation, mixing length
  flowField.getNu(i, j) = 0;
  flowField.getU(i, j) = 0;
  flowField.getLm(i, j) = 0;
}

void NutStencilL::apply(FlowField& flowField, int i, int j, int k) {
  // no vortex viscosity, velociy fluctuation, mixing length
  flowField.getNu(i, j, k) = 0;
  flowField.getU(i, j, k) = 0;
  flowField.getLm(i, j, k) = 0;
}

class limiter1 : public limiter {
 private:
 public:
  limiter1() {}

  virtual ~limiter1() {}

  void limit(const Parameters& parameters, FLOAT& lm, int i, int j) {}

  void limit(const Parameters& parameters, FLOAT& lm, int i, int j, int k) {}
};

class limiter3 : public limiter {
 private:
 public:
  limiter3() {}

  virtual ~limiter3() {}

  void limit(const Parameters& parameters, FLOAT& lm, int i, int j) {
    // only valid for channel!
    FLOAT x = parameters.meshsize->getDx(i, j) / 2 +
              parameters.meshsize->getPosX(i, j);
    FLOAT Rex = x * parameters.flow.Re * parameters.walls.vectorLeft[0];
    lm = std::min(lm, 0.09 * 4.91 * x / sqrt(Rex));
  }

  void limit(const Parameters& parameters, FLOAT& lm, int i, int j, int k) {
    // only valid for channel!
    FLOAT x = parameters.meshsize->getDx(i, j, k) / 2 +
              parameters.meshsize->getPosX(i, j, k);
    FLOAT Rex = x * parameters.flow.Re * parameters.walls.vectorLeft[0];
    lm = std::min(lm, 0.09 * 4.91 * x / sqrt(Rex));
  }
};

class limiter4 : public limiter {
 private:
 public:
  limiter4() {}

  virtual ~limiter4() {}

  void limit(const Parameters& parameters, FLOAT& lm, int i, int j) {
    // only valid for channel!
    FLOAT x = parameters.meshsize->getDx(i, j) / 2 +
              parameters.meshsize->getPosX(i, j);
    FLOAT Rex = x * parameters.flow.Re * parameters.walls.vectorLeft[0];
    lm = std::min(lm, 0.09 * 0.382 * x / pow(Rex, 0.2));
  }

  void limit(const Parameters& parameters, FLOAT& lm, int i, int j, int k) {
    // only valid for channel!
    FLOAT x = parameters.meshsize->getDx(i, j, k) / 2 +
              parameters.meshsize->getPosX(i, j, k);
    FLOAT Rex = x * parameters.flow.Re * parameters.walls.vectorLeft[0];
    lm = std::min(lm, 0.09 * 0.382 * x / pow(Rex, 0.2));
  }
};

NutStencilA::NutStencilA(const Parameters& parameters)
    : FieldStencil<FlowField>(parameters) {
  std::string nulimiter = parameters.simulation.nulimiter;

  if (nulimiter == "1") {
    // no limitation: lm = kappa*h
    l = new limiter1();
  } else if (nulimiter == "2") {
    // TODO: extract the local boundary layer thickness from the laminar
    // reference case
    l = new limiter1();
  } else if (nulimiter == "3") {
    // compute the local boundary layer thickness by assuming a laminar flat
    // plate Blasius boundary layer
    l = new limiter3();
  } else if (nulimiter == "4") {
    // compute the local boundary layer thickness by assuming a turbulent flat
    // boundary layer
    l = new limiter4();
  } else {
    // default: no limitionen
    l = new limiter1();
  }
}

void NutStencilA::apply(FlowField& flowField, int i, int j) {
  // Load local velocities into the center layer of the local array
  loadLocalVelocity2D(flowField, _localVelocity, i, j);
  loadLocalMeshsize2D(_parameters, _localMeshsize, i, j);

  computeNUT2D(i, j, _localVelocity, _localMeshsize, _parameters,
               flowField.getH(i, j), flowField.getNu(i, j),
               flowField.getU(i, j), flowField.getLm(i, j));
}

void NutStencilA::apply(FlowField& flowField, int i, int j, int k) {
  // Load local velocities into the center layer of the local array
  loadLocalVelocity3D(flowField, _localVelocity, i, j, k);
  loadLocalMeshsize3D(_parameters, _localMeshsize, i, j, k);

  computeNUT3D(i, j, k, _localVelocity, _localMeshsize, _parameters,
               flowField.getH(i, j, k), flowField.getNu(i, j, k),
               flowField.getU(i, j, k), flowField.getLm(i, j, k));
}

void NutStencilA::computeNUT2D(int i, int j, const FLOAT* const localVelocity,
                               const FLOAT* const localMeshsize,
                               const Parameters& parameters, const FLOAT& h,
                               FLOAT& nu, FLOAT& flu, FLOAT& lmm) {
  const FLOAT grad = sqrt(2 * computeSijSij2D(localVelocity, localMeshsize));

  // calculate mixing length
  FLOAT lm = 0.41 * h;
  l->limit(parameters, lm, i, j);
  lmm = lm;

  // calculate velocity fluctuation
  flu = lm * grad;

  // calculate vortex viscosity
  nu = lm * lm * grad;
}

void NutStencilA::computeNUT3D(int i, int j, int k,
                               const FLOAT* const localVelocity,
                               const FLOAT* const localMeshsize,
                               const Parameters& parameters, const FLOAT& h,
                               FLOAT& nu, FLOAT& flu, FLOAT& lmm) {
  const FLOAT grad = sqrt(2 * computeSijSij3D(localVelocity, localMeshsize));

  // calculate mixing length
  FLOAT lm = 0.41 * h;
  l->limit(parameters, lm, i, j, k);
  lmm = lm;

  // calculate velocity fluctuation
  flu = lm * grad;

  // calculate vortex viscosity
  nu = lm * lm * grad;
}
}
}
}
