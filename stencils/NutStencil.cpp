#include <algorithm>

#include "NutStencil.h"
#include "StencilFunctions.h"
#include "Definitions.h"

NutStencil::NutStencil(const Parameters& parameters)
    : FieldStencil<FlowFieldTurbA>(parameters) {
  if (parameters.flow.type == "laminar") {
    _fs = new NutStencilL(parameters);
  } else {
    _fs = new NutStencilA(parameters);
  }
}

NutStencil::~NutStencil() { delete _fs; }

void NutStencil::apply(FlowFieldTurbA& flowField, int i, int j) {
  _fs->apply(flowField, i, j);
}

void NutStencil::apply(FlowFieldTurbA& flowField, int i, int j, int k) {
  _fs->apply(flowField, i, j, k);
}

NutStencilL::NutStencilL(const Parameters& parameters)
    : FieldStencil<FlowFieldTurbA>(parameters) {}

void NutStencilL::apply(FlowFieldTurbA& flowField, int i, int j) {
  flowField.getNu(i, j) = 0;
  flowField.getU(i, j) = 0;
  flowField.getLm(i, j) = 0;
}

void NutStencilL::apply(FlowFieldTurbA& flowField, int i, int j, int k) {
  flowField.getNu(i, j, k) = 0;
  flowField.getU(i, j, k) = 0;
  flowField.getLm(i, j, k) = 0;
}

NutStencilA::NutStencilA(const Parameters& parameters)
    : FieldStencil<FlowFieldTurbA>(parameters) {}

void NutStencilA::apply(FlowFieldTurbA& flowField, int i, int j) {
  // Load local velocities into the center layer of the local array
  loadLocalVelocity2D(flowField, _localVelocity, i, j);
  loadLocalMeshsize2D(_parameters, _localMeshsize, i, j);

  computeNUT2D(_localVelocity, _localMeshsize, _parameters,
               flowField.getH(i, j), flowField.getNu(i, j),
               flowField.getU(i, j), flowField.getLm(i, j));
}

void NutStencilA::apply(FlowFieldTurbA& flowField, int i, int j, int k) {
  // Load local velocities into the center layer of the local array
  loadLocalVelocity3D(flowField, _localVelocity, i, j, k);
  loadLocalMeshsize3D(_parameters, _localMeshsize, i, j, k);

  computeNUT3D(_localVelocity, _localMeshsize, _parameters,
               flowField.getH(i, j, k), flowField.getNu(i, j, k),
               flowField.getU(i, j, k), flowField.getLm(i, j, k));
}
