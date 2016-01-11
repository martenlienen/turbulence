#ifndef _NSEOF_FLOWMODELS_TURBULENT_FLOWFIELD_H_
#define _NSEOF_FLOWMODELS_TURBULENT_FLOWFIELD_H_

#include "../../DataStructures.h"
#include "../../Parameters.h"
#include "../../FlowField.h"

namespace nseof {

namespace flowmodels {

namespace turbulent {

/** Flow field for a general turbulence model:
 *    - it is a flow field
 *    - has additional scalar fields
 *        (vortex viscosity, wall distance, velocity fluctuation, mixing length)
 */
class FlowField : public nseof::FlowField {
 public:
  FlowField(const Parameters& parameters);

  // vortex viscosity
  FLOAT& getNu(int i, int j);
  FLOAT& getNu(int i, int j, int k);

  // wall distance
  FLOAT& getH(int i, int j);
  FLOAT& getH(int i, int j, int k);

  // velocity fluctuation
  FLOAT& getU(int i, int j);
  FLOAT& getU(int i, int j, int k);

  // mixing length
  FLOAT& getLm(int i, int j);
  FLOAT& getLm(int i, int j, int k);

 private:
  ScalarField _nu;
  ScalarField _h;
  ScalarField _u;
  ScalarField _lm;
};
}
}
}

#endif
