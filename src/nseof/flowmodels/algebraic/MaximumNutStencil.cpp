#include "MaximumNutStencil.h"

#include <algorithm>

namespace nseof {

namespace flowmodels {

namespace algebraic {

MaximumNutStencil::MaximumNutStencil(const Parameters& parameters)
    : FieldStencil<FlowField>(parameters) {
  this->maximum = 0;
}

MaximumNutStencil::~MaximumNutStencil() {}

void MaximumNutStencil::apply(FlowField& flowField, int i, int j) {
  this->maximum = std::max(maximum, flowField.getNu(i, j));
}

void MaximumNutStencil::apply(FlowField& flowField, int i, int j, int k) {
  this->maximum = std::max(maximum, flowField.getNu(i, j, k));
}

void MaximumNutStencil::reset() {
  this->maximum = 0;
}
}
}
}
