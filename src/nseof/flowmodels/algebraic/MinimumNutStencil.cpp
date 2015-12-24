#include "MinimumNutStencil.h"

#include <algorithm>

namespace nseof {

namespace flowmodels {

namespace algebraic {

MinimumNutStencil::MinimumNutStencil(const Parameters& parameters)
    : FieldStencil<FlowField>(parameters) {
  this->minimum = 0;
}

MinimumNutStencil::~MinimumNutStencil() {}

void MinimumNutStencil::apply(FlowField& flowField, int i, int j) {
  this->minimum = std::max(minimum, flowField.getNu(i, j));
}

void MinimumNutStencil::apply(FlowField& flowField, int i, int j, int k) {
  this->minimum = std::max(minimum, flowField.getNu(i, j, k));
}

void MinimumNutStencil::reset() {
  this->minimum = 0;
}
}
}
}
