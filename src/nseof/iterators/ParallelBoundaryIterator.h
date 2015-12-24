#ifndef _PARALLEL_BOUNDARY_ITERATOR_H_
#define _PARALLEL_BOUNDARY_ITERATOR_H_

#include <vector>

#include "../Iterators.h"

#include "RangeIterator.h"

namespace nseof {

template <typename FF, typename T>
class ParallelBoundaryIterator : Iterator<FF> {
 public:
  RangeIterator<FF, T> first;
  RangeIterator<FF, T> second;

  ParallelBoundaryIterator(
      FF& flowField, const Parameters& parameters, Point lowStart, Point lowEnd,
      Point highStart, Point highEnd,
      std::function<void(FF& flowField, int, int, int, T&)> apply)
      : Iterator<FF>(flowField, parameters),
        first(flowField, parameters, lowStart, lowEnd, apply),
        second(flowField, parameters, highStart, highEnd, apply){};

  void iterate();

  /**
   * Convert an offset to the right-top-back corner to a point relative to the
   * origin.
   */
  static Point highOffsetToPoint(FF& flowField, Point offset) {
    return {flowField.getCellsX() + offset.x - 1,
            flowField.getCellsY() + offset.y - 1,
            flowField.getCellsZ() + offset.z - 1};
  }
};

template <typename FF, typename T>
void ParallelBoundaryIterator<FF, T>::iterate() {
  this->first.iterate();
  this->second.iterate();
}
}

#endif  // _PARALLEL_BOUNDARY_ITERATOR_H_
