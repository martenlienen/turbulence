#include "HStencil.h"

namespace nseof {

namespace flowmodels {

namespace algebraic {

HStencil::HStencil(const Parameters& parameters)
    : FieldStencil<FlowField>(parameters), wdm(_parameters) {
  wdm.init();
}

void HStencil::apply(FlowField& flowField, int i, int j) {
  const int obstacle = flowField.getFlags().getValue(i, j);

  if ((obstacle & OBSTACLE_SELF) == 0) {
    flowField.getH(i, j) = wdm.query(i, j);
  } else {
    flowField.getH(i, j) = 0.0;
  }
}

void HStencil::apply(FlowField& flowField, int i, int j, int k) {
  const int obstacle = flowField.getFlags().getValue(i, j, k);

  if ((obstacle & OBSTACLE_SELF) == 0) {
    flowField.getH(i, j, k) = wdm.query(i, j, k);
  } else {
    flowField.getH(i, j, k) = 0.0;
  }
}
}
}
}
