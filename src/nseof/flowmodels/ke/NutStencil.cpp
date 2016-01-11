#include <algorithm>

#include "../../stencils/StencilFunctions.h"

#include "NutStencil.h"

namespace nseof {

namespace flowmodels {

namespace ke {

NutStencil::NutStencil(const Parameters& parameters)
    : FieldStencil<FlowField>(parameters) {
  // do you want to calculate vortex viscosity?
  if (parameters.flow.type == "laminar") {
    // no
    _fs = new NutStencilL(parameters);
  } else {
    // yes
    _fs = new NutStencilA(parameters, parameters.kEpsilon.cmu);
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
}

void NutStencilL::apply(FlowField& flowField, int i, int j, int k) {
  // no vortex viscosity, velociy fluctuation, mixing length
  flowField.getNu(i, j, k) = 0;
  flowField.getU(i, j, k) = 0;
}

NutStencilA::NutStencilA(const Parameters& parameters, FLOAT cmu)
    : FieldStencil<FlowField>(parameters), cmu(cmu) {}

void NutStencilA::apply(FlowField& flowField, int i, int j) {
  apply(flowField, i, j, 0);
}

void NutStencilA::apply(FlowField& flowField, int i, int j, int k) {
  FLOAT tke = fabs(flowField.getTke(i, j, k));
  FLOAT e = fabs(flowField.getEpsilon(i, j, k));

  FLOAT nut = cmu * flowField.getFmu(i, j, k) * tke * tke / (max(e, 1e-8));
  nut = max(nut, 1E-4 / _parameters.flow.Re);
  nut = min(nut, 1E+2 / _parameters.flow.Re);

  flowField.getNu(i, j, k) = nut;
  flowField.getLm(i, j, k) = cmu * flowField.getFmu(i, j, k) * pow(tke, 1.5) *
                             sqrt(1.5) / (max(e, 1e-8));
  flowField.getU(i, j, k) = sqrt(2 / 3 * tke);
}
}
}
}
