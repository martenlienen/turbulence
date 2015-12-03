#ifndef _BOTTOM_TOP_BOUNDARY_ITERATOR_H_
#define _BOTTOM_TOP_BOUNDARY_ITERATOR_H_

#include "ParallelBoundaryIterator.h"

template <typename FF, typename T>
class BottomTopBoundaryIterator : public ParallelBoundaryIterator<FF, T> {
 public:
  BottomTopBoundaryIterator(
      FF& flowField, const Parameters& parameters, Point lowOffset,
      Point highOffset,
      std::function<void(FF& flowField, int, int, int, T&)> apply)
      : BottomTopBoundaryIterator(
            flowField, parameters, lowOffset,
            ParallelBoundaryIterator<FF, T>::highOffsetToPoint(flowField,
                                                               highOffset),
            apply, false) {}

 private:
  BottomTopBoundaryIterator(
      FF& flowField, const Parameters& parameters, Point low, Point high,
      std::function<void(FF& flowField, int, int, int, T&)> apply, bool x)
      : ParallelBoundaryIterator<FF, T>(
            flowField, parameters, low, Point(high.x, low.y, high.z),
            Point(low.x, high.y, low.z), high, apply){};
};

#endif  // _BOTTOM_TOP_BOUNDARY_ITERATOR_H_
