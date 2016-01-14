#include <ctime>
#include <fstream>
#include <functional>
#include <iostream>

#include "../Definitions.h"

#include "GeometryManager.h"

namespace nseof {

namespace geometry {

GeometryManager::GeometryManager(const Parameters& parameters)
    : _parameters(parameters), x(0), y(0), z(0) {
  if (parameters.geometry.obstacle != "") {
    std::string file(parameters.geometry.obstacle);

    std::ifstream infile(file);

    int a, b, c;

    while (infile >> a >> b >> c) {
      x.push_back(a);
      y.push_back(b);
      z.push_back(c);
    }
  }
}

GeometryManager::~GeometryManager() {}

void GeometryManager::init(FlowField& flowField) {
  IntScalarField& f = flowField.getFlags();

  for (unsigned int i = 0; i < x.size(); i++) {
    int xx = x[i] + 2;
    int yy = y[i] + 2;
    int zz = _parameters.geometry.dim == 2 ? 0 : z[i] + 2;

    f.getValue(xx, yy, zz) |= OBSTACLE_SELF;

    f.getValue(xx - 1, yy + 0, zz + 0) |= OBSTACLE_RIGHT;
    f.getValue(xx - 1, yy + 0, zz + 0) |= OBSTACLE_LEFT;
    f.getValue(xx + 0, yy - 1, zz + 0) |= OBSTACLE_TOP;
    f.getValue(xx + 0, yy + 1, zz + 0) |= OBSTACLE_BOTTOM;

    if (_parameters.geometry.dim == 3) {
      f.getValue(xx + 0, yy + 0, zz - 1) |= OBSTACLE_BACK;
      f.getValue(xx + 0, yy + 0, zz + 1) |= OBSTACLE_FRONT;
    }
  }
}
}
}
