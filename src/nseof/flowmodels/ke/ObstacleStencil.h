#ifndef _NSEOF_FLOWMODELS_KE_OBSTACLE_STENCIL_H_
#define _NSEOF_FLOWMODELS_KE_OBSTACLE_STENCIL_H_

#include "../../Parameters.h"
#include "../../Stencil.h"

#include "FlowField.h"
#include "BFInputStencils.h"

namespace nseof {

namespace flowmodels {

namespace ke {

/**
 * Compute all velocities on obstacle cells. This has been taken out of
 * Velocity stencil to circumvent any race conditions
 */
class ObstacleStencil : public FieldStencil<FlowField> {
 private:
  WallStencil _wke;

 public:
  /* Constructor
   * @param parameters Parameters of the problem
   */
  ObstacleStencil(const Parameters& parameters);

  /** Apply the stencil in 2D
   * @param flowField Flow field information
   * @param i Position in the X direction
   * @param j Position in the Y direction
   */
  void apply(FlowField& flowField, int i, int j);

  /** Apply the stencil in 3D
   * @param flowField Flow field information
   * @param i Position in the X direction
   * @param j Position in the Y direction
   * @param k Position in the Z direction
   */
  void apply(FlowField& flowField, int i, int j, int k);
};
}
}
}

#endif
