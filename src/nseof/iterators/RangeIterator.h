#ifndef RANGE_ITERATOR_H_
#define RANGE_ITERATOR_H_

#include <functional>
#include <vector>

#include "../Iterators.h"
#include "../Point.h"

namespace nseof {

template <typename FF, typename T>
class RangeIterator : public Iterator<FF> {
 public:
  RangeIterator(FF& flowField, const Parameters& parameters, Point low,
                Point high,
                std::function<void(FF& flowField, int, int, int, T&)> apply)
      : Iterator<FF>(flowField, parameters),
        flowField(flowField),
        low(low),
        high(high),
        apply(apply),
        data((high - low).size()) {}

  void iterate();

 private:
  FF& flowField;
  Point low;
  Point high;
  std::function<void(FF& flowField, int, int, int, T&)> apply;

 public:
  std::vector<T> data;
};

template <typename FF, typename T>
void RangeIterator<FF, T>::iterate() {
  Point& l = this->low;
  Point& h = this->high;
  Point diff = h - l;
  diff.x++;
  diff.y++;
  diff.z++;

  if (this->_parameters.geometry.dim == 2) {
    for (int i = l.x; i <= h.x; i++) {
      for (int j = l.y; j <= h.y; j++) {
        int m = (j - l.y) + (i - l.x) * diff.y;

        this->apply(this->flowField, i, j, 0, this->data[m]);
      }
    }
  } else {
    for (int i = l.x; i <= h.x; i++) {
      for (int j = l.y; j <= h.y; j++) {
        for (int k = l.z; k <= h.z; k++) {
          int m =
              (k - l.z) + (j - l.y) * diff.z + (i - l.x) * (diff.y * diff.z);

          this->apply(this->flowField, i, j, k, this->data[m]);
        }
      }
    }
  }
}
}

#endif  // RANGE_ITERATOR_H_
