#ifndef _NSEOF_FLOWMODELS_KE_RHS_STENCIL_H_
#define _NSEOF_FLOWMODELS_KE_RHS_STENCIL_H_

#include "../../Definitions.h"
#include "../../Parameters.h"
#include "../../Stencil.h"

#include "FlowField.h"

namespace nseof {

namespace flowmodels {

namespace ke {

class RHSStencil : public FieldStencil<FlowField> {
 private:
  // A local velocity variable that will be used to approximate derivatives.
  // Size matches 3D
  // case, but can be used for 2D as well.
  FLOAT _localVelocity[27 * 3];
  // local meshsize
  FLOAT _localMeshsize[27 * 3];
  // pm: 2015.11.07
  // local nu
  FLOAT _localNu[27 * 3];
  // local tke
  FLOAT _localTKE[27 * 3];
  // local epsilon
  FLOAT _localEpsilon[27 * 3];
  // local fmu*nuT
  FLOAT _localFmuNut[27 * 3];

 public:
  RHSStencil(const Parameters& parameters);

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
