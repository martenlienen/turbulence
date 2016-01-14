#ifndef _NSEOF_GEOMETRY_GEOMETRYMANAGER_H_
#define _NSEOF_GEOMETRY_GEOMETRYMANAGER_H_

#include <ANN/ANN.h>
#include <vector>

#include "../Parameters.h"
#include "../FlowField.h"

namespace nseof {

namespace geometry {

class GeometryManager {
 private:
  // general field
  const Parameters& _parameters;

 public:
  std::vector<double> x;
  std::vector<double> y;
  std::vector<double> z;

  GeometryManager(const Parameters& parameters);

  virtual ~GeometryManager();

  void init(FlowField& flowField);
};
}
}

#endif  // _NSEOF_GEOMETRY_GEOMETRYMANAGER_H_
