#ifndef _NSEOF_FLOWMODELS_ALGEBRAIC_FGHSTENCIL_H_
#define _NSEOF_FLOWMODELS_ALGEBRAIC_FGHSTENCIL_H_

#include "../../Stencil.h"
#include "../../Parameters.h"

#include "FlowField.h"

namespace nseof {

namespace flowmodels {

namespace algebraic {

class FGHStencil : public FieldStencil<FlowField> {
 private:
  // A local velocity variable that will be used to approximate derivatives.
  // Size matches 3D case, but can be used for 2D as well.
  FLOAT _localVelocity[27 * 3];
  // local meshsize
  FLOAT _localMeshsize[27 * 3];
  // local nu
  FLOAT _localNu[27 * 3];

 public:
  FGHStencil(const Parameters& parameters);

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
