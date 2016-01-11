#include "../../stencils/StencilFunctions.h"

#include "KEStencil.h"

namespace nseof {

namespace flowmodels {

namespace ke {

KEStencil::KEStencil(const Parameters& parameters)
    : FieldStencil<FlowField>(parameters) {}

void KEStencil::apply(FlowField& flowField, int i, int j) {
  // calculate new value from old value and rhs
  flowField.getTke(i, j) += _parameters.timestep.dt * flowField.getRHSTke(i, j);
  flowField.getTke(i, j) *= flowField.getTke(i, j) < 0.0 ? 0.0 : 1.0;

  flowField.getEpsilon(i, j) +=
      _parameters.timestep.dt * flowField.getRHSEpsilon(i, j);
  flowField.getEpsilon(i, j) *= flowField.getEpsilon(i, j) < 0.0 ? 0.0 : 1.0;
}

void KEStencil::apply(FlowField& flowField, int i, int j, int k) {
  // calculate new value from old value and rhs
  flowField.getTke(i, j, k) +=
      _parameters.timestep.dt * flowField.getRHSTke(i, j, k);
  flowField.getTke(i, j, k) *= flowField.getTke(i, j, k) < 0.0 ? 0.0 : 1.0;

  flowField.getEpsilon(i, j, k) +=
      _parameters.timestep.dt * flowField.getRHSEpsilon(i, j, k);
  flowField.getEpsilon(i, j, k) *=
      flowField.getEpsilon(i, j, k) < 0.0 ? 0.0 : 1.0;
}
}
}
}
