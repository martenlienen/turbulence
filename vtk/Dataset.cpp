#include "Dataset.h"

#include <ios>
#include <sstream>

namespace vtk {

Dataset::Dataset(Point cells, Point points, std::vector<Point> coordinates)
    : cells(cells), points(points), coordinates(coordinates) {}

int Dataset::numCells() {
  Point c = this->cells;

  return c.x * c.y * c.z;
}

std::string Dataset::str() {
  Point p = this->points;

  std::ostringstream stream;

  stream << std::fixed;

  stream << "DATASET STRUCTURED_GRID\n";
  stream << "DIMENSIONS " << (int)p.x << " " << (int)p.y << " " << (int)p.z
         << "\n";
  stream << "POINTS " << (int)(p.x * p.y * p.z) << " float\n";

  for (auto& c : this->coordinates) {
    stream << c.x << " " << c.y << " " << c.z << "\n";
  }

  return stream.str();
}
}
