#ifndef _NSEOF_FLOWMODELS_KE_KE_STENCIL_F_H_
#define _NSEOF_FLOWMODELS_KE_KE_STENCIL_F_H_

#include <functional>

#include "../../Definitions.h"
#include "../../Parameters.h"
#include "../../Stencil.h"
#include "../../walldistance/WallDistanceManager.h"

#include "FlowField.h"

namespace nseof {

namespace flowmodels {

namespace ke {

class KEStencilF : public FieldStencil<FlowField> {
 private:
  FLOAT err = 1e-7;

  struct woparameters {
    FLOAT fmu = 1.0;
    FLOAT f1 = 1.0;
    FLOAT f2 = 1.0;
    FLOAT f3 = 1.0;
    FLOAT D = 0.0;
    FLOAT E = 0.0;
  } wop;

  struct wiparameters {
    FLOAT Re = 1.0;
    FLOAT tke = 1.0;
    FLOAT tkes = 1.0;
    FLOAT epsilon = 1.0;
    FLOAT epsilons = 1.0;
    FLOAT delta = 0.0;
    FLOAT Rt = 0.0;
    FLOAT Rd = 0.0;
    FLOAT yplus = 0.0;
  } wip;

  std::function<void()> calc;

 public:
  FLOAT cmu;

  KEStencilF(const Parameters& parameters);

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
