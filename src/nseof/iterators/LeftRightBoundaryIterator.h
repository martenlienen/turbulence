#ifndef _LEFT_RIGHT_BOUNDARY_ITERATOR_H_
#define _LEFT_RIGHT_BOUNDARY_ITERATOR_H_

#include "ParallelBoundaryIterator.h"

namespace nseof {

template <typename FF, typename T>
class LeftRightBoundaryIterator : public ParallelBoundaryIterator<FF, T> {
 public:
  LeftRightBoundaryIterator(
      FF& flowField, const Parameters& parameters, Point lowOffset,
      Point highOffset,
      std::function<void(FF& flowField, int, int, int, T&)> apply,
      int leftWidth = 1, int rightWidth = 1)
      : LeftRightBoundaryIterator(
            flowField, parameters, lowOffset,
            ParallelBoundaryIterator<FF, T>::highOffsetToPoint(flowField,
                                                               highOffset),
            apply, leftWidth, rightWidth, false) {}

 private:
  LeftRightBoundaryIterator(
      FF& flowField, const Parameters& parameters, Point low, Point high,
      std::function<void(FF& flowField, int, int, int, T&)> apply,
      int leftWidth, int rightWidth, bool x)
      : ParallelBoundaryIterator<FF, T>(
            flowField, parameters, low,
            Point(low.x + (leftWidth - 1), high.y, high.z),
            Point(high.x - (rightWidth - 1), low.y, low.z), high, apply){};
};
}

#endif  // _LEFT_RIGHT_BOUNDARY_ITERATOR_H_
