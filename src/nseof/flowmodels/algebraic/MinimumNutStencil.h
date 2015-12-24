#ifndef _NSEOF_FLOWMODELS_ALGEBRAIC_MINIMUMNUTSTENCIL_H_
#define _NSEOF_FLOWMODELS_ALGEBRAIC_MINIMUMNUTSTENCIL_H_

#include "../../Definitions.h"
#include "../../Parameters.h"
#include "../../Stencil.h"

#include "FlowField.h"

namespace nseof {

namespace flowmodels {

namespace algebraic {

/**
 * When iterated with, finds and stores the minimum nu_t value
 */
class MinimumNutStencil : public FieldStencil<FlowField> {
 public:
  MinimumNutStencil(const Parameters& parameters);

  ~MinimumNutStencil();

  /** 2D operation for one position
   *
   * @param flowField State of the flow field
   * @param i Position in the x direction
   * @param j Position in the y direction
   */
  void apply(FlowField& flowField, int i, int j);

  /** 3D operation for one position
   *
   * @param flowField State of the flow field
   * @param i Position in the x direction
   * @param j Position in the y direction
   * @param k Position in the z direction
   */
  void apply(FlowField& flowField, int i, int j, int k);

  double getMinimum() { return minimum; }

  void reset();

 private:
  double minimum;
};
}
}
}

#endif  // _NSEOF_FLOWMODELS_ALGEBRAIC_MINIMUMNUTSTENCIL_H_
