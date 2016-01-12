#ifndef _NSEOF_FLOWMODELS_KE_KE_STENCIL_H_
#define _NSEOF_FLOWMODELS_KE_KE_STENCIL_H_

#include "../../Stencil.h"
#include "../../Parameters.h"
#include "../../walldistance/WallDistanceManager.h"

#include "FlowField.h"

namespace nseof {

namespace flowmodels {

namespace ke {

class KEStencil : public FieldStencil<FlowField> {
 public:
  KEStencil(const Parameters& parameters);

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
};
}
}
}

#endif
