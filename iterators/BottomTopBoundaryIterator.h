#ifndef _BOTTOM_TOP_BOUNDARY_ITERATOR_H_
#define _BOTTOM_TOP_BOUNDARY_ITERATOR_H_

#include "ParallelBoundaryIterator.h"

template <typename FF, typename T>
class BottomTopBoundaryIterator : public ParallelBoundaryIterator<FF, T> {
 public:
  BottomTopBoundaryIterator(
      FF& flowField, const Parameters& parameters, Point lowOffset,
      Point highOffset,
      std::function<void(FF& flowField, int, int, int, T&)> apply,
      int bottomWidth = 1, int topWidth = 1)
      : BottomTopBoundaryIterator(
            flowField, parameters, lowOffset,
            ParallelBoundaryIterator<FF, T>::highOffsetToPoint(flowField,
                                                               highOffset),
            apply, bottomWidth, topWidth, false) {}

 private:
  BottomTopBoundaryIterator(
      FF& flowField, const Parameters& parameters, Point low, Point high,
      std::function<void(FF& flowField, int, int, int, T&)> apply,
      int bottomWidth, int topWidth, bool x)
      : ParallelBoundaryIterator<FF, T>(
            flowField, parameters, low,
            Point(high.x, low.y + (bottomWidth - 1), high.z),
            Point(low.x, high.y - (bottomWidth - 1), low.z), high, apply){};
};

#endif  // _BOTTOM_TOP_BOUNDARY_ITERATOR_H_
