#ifndef _NSEOF_FLOWMODELS_ALGEBRAIC_MAXIMUMNUTSTENCIL_H_
#define _NSEOF_FLOWMODELS_ALGEBRAIC_MAXIMUMNUTSTENCIL_H_

#include "../../Definitions.h"
#include "../../Parameters.h"
#include "../../Stencil.h"

#include "FlowField.h"

namespace nseof {

namespace flowmodels {

namespace algebraic {

/**
 * When iterated with, finds and stores the maximum nu_t value
 */
class MaximumNutStencil : public FieldStencil<FlowField> {
 public:
  MaximumNutStencil(const Parameters& parameters);

  ~MaximumNutStencil();

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

  double getMaximum() { return maximum; }

  void reset();

 private:
  double maximum;
};
}
}
}

#endif  // _NSEOF_FLOWMODELS_ALGEBRAIC_MAXIMUMNUTSTENCIL_H_
