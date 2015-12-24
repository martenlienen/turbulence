#ifndef _NSEOF_FLOWMODELS_ALGEBRAIC_FLOWFIELD_H_
#define _NSEOF_FLOWMODELS_ALGEBRAIC_FLOWFIELD_H_

#include "../../DataStructures.h"
#include "../../Parameters.h"
#include "../../FlowField.h"

namespace nseof {

namespace flowmodels {

namespace algebraic {

/** Flow field for algebraic turbulent model:
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

// Load the local viscosity (nu+nut) cube with relevant viscosity of the 2D
// plane
void loadLocalNu2D(const Parameters& parameters, FlowField& flowField,
                   FLOAT* const localNu, int i, int j);

// Load the local viscosity (nu+nut) cube with surrounding values
void loadLocalNu3D(const Parameters& parameters, FlowField& flowField,
                   FLOAT* const localNu, int i, int j, int k);
}
}
}

#endif
