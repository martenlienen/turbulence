#ifndef _POINT_H_
#define _POINT_H_

#include <initializer_list>
#include <vector>

namespace nseof {

struct Point {
  int x;
  int y;
  int z;

  Point(int x, int y, int z) : x(x), y(y), z(z) {}
  Point(std::initializer_list<int> list) {
    std::vector<int> coords(list);
    int size = list.size();

    if (size > 0) {
      this->x = coords[0];
    }

    if (size > 1) {
      this->y = coords[1];
    }

    if (size > 2) {
      this->z = coords[2];
    }
  };

  Point operator-(const Point& other) {
    return Point(this->x - other.x, this->y - other.y, this->z - other.z);
  }

  // TODO comment
  int size() { return (this->x + 1) * (this->y + 1) * (this->z + 1); }
};
}

#endif  // _POINT_H_
