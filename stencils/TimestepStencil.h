#ifndef _TIMESTEP_STENCIL_P_H_
#define _TIMESTEP_STENCIL_P_H_

#include "../Definitions.h"
#include "../Parameters.h"
#include "../Stencil.h"
#include "../FlowFieldTurbA.h"
#include <string>
#include <fstream>
#include <sstream>
#include <iostream>
using namespace std;

/** TODO WS1: Stencil for writting VTK files
 *
 * When iterated with, creates a VTK file.
 */
class TimestepStencil : public FieldStencil<FlowFieldTurbA> {
 public:
  /** Constructor
   *
   * @param prefix String with the prefix of the name of the VTK files
   */
  TimestepStencil(const Parameters& parameters);

  ~TimestepStencil();

  /** 2D operation for one position
   *
   * @param flowField State of the flow field
   * @param i Position in the x direction
   * @param j Position in the y direction
   */
  void apply(FlowFieldTurbA& flowField, int i, int j);

  /** 3D operation for one position
   *
   * @param flowField State of the flow field
   * @param i Position in the x direction
   * @param j Position in the y direction
   * @param k Position in the z direction
   */
  void apply(FlowFieldTurbA& flowField, int i, int j, int k);

  double getMinimum() { return minimum; }

 private:
  double minimum;
};

#endif
