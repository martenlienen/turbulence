#ifndef _NSEOF_FLOWMODELS_ALGEBRAIC_HSTENCIL_H_
#define _NSEOF_FLOWMODELS_ALGEBRAIC_HSTENCIL_H_

#include "../../Stencil.h"
#include "../../Parameters.h"
#include "../../walldistance/WallDistanceManager.h"

#include "FlowField.h"

namespace nseof {

namespace flowmodels {

namespace algebraic {

class HStencil : public FieldStencil<FlowField> {
 public:
  HStencil(const Parameters& parameters);

  /** Apply the stencil in 2D
   *
   * Performs the operation of the stencil in a single position given by the
   * indexes.
   * @param flowField State of the flow
   * @param i Index in the x direction
   * @param j Index in the y direction
   */
  void apply(FlowField& flowField, int i, int j);

  /** Apply the stencil in 3D
   *
   * @param flowField State of the flow
   * @param i Index in the x direction
   * @param j Index in the y direction
   * @param k Index in the z direction
   */
  void apply(FlowField& flowField, int i, int j, int k);

 private:
  WallDistanceManager wdm;
};
}
}
}

#endif
