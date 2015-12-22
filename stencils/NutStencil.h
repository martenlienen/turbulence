#ifndef _STENCIL_NUT_H_
#define _STENCIL_NUT_H_

#include "../FlowFieldTurbA.h"
#include "../Stencil.h"
#include "../Parameters.h"

class NutStencil : public FieldStencil<FlowFieldTurbA> {
 private:
  FLOAT _localVelocity[27 * 3];
  FLOAT _localMeshsize[27 * 3];

 public:
  NutStencil(const Parameters& parameters);

  /** Apply the stencil in 2D
   *
   * Performs the operation of the stencil in a single position given by the
   * indexes.
   * @param flowField State of the flow
   * @param i Index in the x direction
   * @param j Index in the y direction
   */
  void apply(FlowFieldTurbA& flowField, int i, int j);

  /** Apply the stencil in 3D
   *
   * @param flowField State of the flow
   * @param i Index in the x direction
   * @param j Index in the y direction
   * @param k Index in the z direction
   */
  void apply(FlowFieldTurbA& flowField, int i, int j, int k);
};

#endif
