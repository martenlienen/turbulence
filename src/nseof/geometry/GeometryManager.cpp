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
  for (unsigned long int i = 0; i < x.size(); i++) {
    flowField.getFlags().getValue(x[i] + 2, y[i] + 2, z[i] + 2) = OBSTACLE_SELF;
  }

  for (unsigned long int i = 0; i < x.size(); i++) {
    int xx = x[i] + 2;
    int yy = y[i] + 2;
    int zz = _parameters.geometry.dim == 2 ? 0 : z[i] + 2;

    IntScalarField& f = flowField.getFlags();

    // left cell
    write(f, xx, yy, zz, -1, +0, +0, OBSTACLE_LEFT);
    // right cell
    write(f, xx, yy, zz, +1, +0, +0, OBSTACLE_RIGHT);
    // bottom cell
    write(f, xx, yy, zz, +0, -1, +0, OBSTACLE_BOTTOM);
    // top cell
    write(f, xx, yy, zz, +0, +1, +0, OBSTACLE_TOP);

    if (_parameters.geometry.dim == 3) {
      // bottom cell
      write(f, xx, yy, zz, +0, +0, -1, OBSTACLE_FRONT);
      // top cell
      write(f, xx, yy, zz, +0, +0, +1, OBSTACLE_BACK);
    }
  }
}

void GeometryManager::write(IntScalarField& f, int xx, int yy, int zz, int a,
                            int b, int c, int d) {
  if ((f.getValue(xx + a, yy + b, zz + c) & OBSTACLE_SELF) == 1) {
    f.getValue(xx + 0, yy + 0, zz + 0) += d;
  }
  f.getValue(xx - a, yy - b, zz - c) += d;
}
}
}
