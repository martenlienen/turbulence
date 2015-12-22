#ifndef _STENCIL_H_H_
#define _STENCIL_H_H_

#include "../FlowFieldTurbA.h"
#include "../Stencil.h"
#include "../Parameters.h"
#include "../walldistance/WallDistanceManager.h"

class HStencil : public FieldStencil<FlowFieldTurbA> {
 public:
  HStencil(const Parameters& parameters);

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

 private:
  WallDistanceManager wdm;
};

#endif
