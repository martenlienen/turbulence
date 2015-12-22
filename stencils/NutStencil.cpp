#include <algorithm>

#include "NutStencil.h"
#include "StencilFunctions.h"
#include "Definitions.h"

NutStencil::NutStencil(const Parameters& parameters)
    : FieldStencil<FlowFieldTurbA>(parameters) {}

void NutStencil::apply(FlowFieldTurbA& flowField, int i, int j) {
  // Load local velocities into the center layer of the local array
  loadLocalVelocity2D(flowField, _localVelocity, i, j);
  loadLocalMeshsize2D(_parameters, _localMeshsize, i, j);

  flowField.getNu(i, j) =
      std::min(0.5 * 0.1, computeNUT2D(_localVelocity, _localMeshsize,
                                       _parameters, flowField.getH(i, j)));

  //    flowField.getNu(i, j) = computeNUT2D(_localVelocity, _localMeshsize,
  //    _parameters, flowField.getH(i, j));
}

void NutStencil::apply(FlowFieldTurbA& flowField, int i, int j, int k) {}
