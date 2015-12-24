#include "MinimumNutStencil.h"

#include <algorithm>

namespace nseof {

MinimumNutStencil::MinimumNutStencil(const Parameters& parameters)
    : FieldStencil<FlowFieldTurbA>(parameters) {
  this->minimum = 0;
}

MinimumNutStencil::~MinimumNutStencil() {}

void MinimumNutStencil::apply(FlowFieldTurbA& flowField, int i, int j) {
  this->minimum = std::max(minimum, flowField.getNu(i, j));
}

void MinimumNutStencil::apply(FlowFieldTurbA& flowField, int i, int j, int k) {
  this->minimum = std::max(minimum, flowField.getNu(i, j, k));
}

void MinimumNutStencil::reset() {
  this->minimum = 0;
}
}
