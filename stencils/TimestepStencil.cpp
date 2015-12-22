#include "TimestepStencil.h"
#include <fstream>
#include <algorithm>

TimestepStencil::TimestepStencil(const Parameters& parameters)
    : FieldStencil<FlowFieldTurbA>(parameters) {
  minimum = 10000000000000000;
}

TimestepStencil::~TimestepStencil() {}

void TimestepStencil::apply(FlowFieldTurbA& flowField, int i, int j) {
  double nu = std::min(flowField.getNu(i, j), 1E-6);

  double dx = _parameters.meshsize->getDx(i, j);
  double dy = _parameters.meshsize->getDy(i, j);

  double dt = (_parameters.flow.Re + 1 / nu) / (1 / dx / dx + 1 / dy / dy);
  double dtt = 1e-4;
  //    double dtt = 1e-3;
  minimum = std::min(std::min(minimum, dt), dtt);
  //    minimum = std::min(minimum,dt);
}

void TimestepStencil::apply(FlowFieldTurbA& flowField, int i, int j, int k) {
  double nu = std::min(flowField.getNu(i, j), 1E-6);

  double dx = _parameters.meshsize->getDx(i, j, k);
  double dy = _parameters.meshsize->getDy(i, j, k);
  double dz = _parameters.meshsize->getDy(i, j, k);

  double dt = (_parameters.flow.Re + 1 / nu) /
              (1 / dx / dx + 1 / dy / dy + 1 / dz / dz);
  double dtt = 1e-4;
  minimum = std::min(std::min(minimum, dt), dtt);
}
