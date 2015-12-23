#include "HStencil.h"
#include "StencilFunctions.h"
#include "Definitions.h"
#include <algorithm>
#include <cmath>

HStencil::HStencil(const Parameters& parameters)
    : FieldStencil<FlowFieldTurbA>(parameters) {}

void HStencil::apply(FlowFieldTurbA& flowField, int i, int j) {
  FLOAT temp = _parameters.meshsize->getPosY(i, j) +
               _parameters.meshsize->getDy(i, j) * 0.5;
  flowField.getH(i, j) = std::min(_parameters.geometry.lengthY - temp, temp);
}

void HStencil::apply(FlowFieldTurbA& flowField, int i, int j, int k) {
  FLOAT tempy = _parameters.meshsize->getPosY(i, j, k) +
                _parameters.meshsize->getDy(i, j, k) * 0.5;
  FLOAT dy = std::min(_parameters.geometry.lengthY - tempy, tempy);

  FLOAT tempz = _parameters.meshsize->getPosZ(i, j, k) +
                _parameters.meshsize->getDz(i, j, k) * 0.5;
  FLOAT dz = std::min(_parameters.geometry.lengthY - tempz, tempz);

  flowField.getH(i, j, k) = std::min(dy, dz);
}
