#ifndef _VTK_DATASET_H
#define _VTK_DATASET_H

#include <string>
#include <vector>

namespace vtk {

struct Point {
 public:
  Point() : x(0), y(0), z(0) {}
  Point(double x, double y, double z) : x(x), y(y), z(z) {}

  double x;
  double y;
  double z;
};

class Dataset {
 public:
  Dataset(Point cells, Point points, std::vector<Point> coordinates);

  int numCells();

  std::string str();

 private:
  Point cells;
  Point points;

  std::vector<Point> coordinates;
};
}

#endif  // _VTK_DATASET_H
