#include "../../stencils/StencilFunctions.h"

#include "InitStencil.h"

namespace nseof {

namespace flowmodels {

namespace ke {

InitStencil::InitStencil(const Parameters& parameters)
    : FieldStencil<FlowField>(parameters) {}

void InitStencil::apply(FlowField& flowField, int i, int j) {
  apply(flowField, i, j, 0);
}

void InitStencil::apply(FlowField& flowField, int i, int j, int k) {
  // no tke and epsilon at the beginning
  flowField.getTke(i, j, k) = 0;
  flowField.getEpsilon(i, j, k) = 0;
}
}
}
}
