#ifndef _NSEOF_FLOWMODELS_KE_BF_INPUT_STENCILS_H_
#define _NSEOF_FLOWMODELS_KE_BF_INPUT_STENCILS_H_

#include "../../Stencil.h"
#include "../../Parameters.h"

#include "FlowField.h"

namespace nseof {

namespace flowmodels {

namespace ke {

/**
 * A stencil to set the input velocity in channel flows. Only implements the
 * applyLeftWall(...) methods.
 */
class BFInputStencil : public BoundaryStencil<FlowField> {
 public:
  BFInputStencil(const Parameters& parameters);

  void applyLeftWall(FlowField& flowField, int i, int j);
  void applyRightWall(FlowField& flowField, int i, int j);
  void applyBottomWall(FlowField& flowField, int i, int j);
  void applyTopWall(FlowField& flowField, int i, int j);

  void applyLeftWall(FlowField& flowField, int i, int j, int k);
  void applyRightWall(FlowField& flowField, int i, int j, int k);
  void applyBottomWall(FlowField& flowField, int i, int j, int k);
  void applyTopWall(FlowField& flowField, int i, int j, int k);
  void applyFrontWall(FlowField& flowField, int i, int j, int k);
  void applyBackWall(FlowField& flowField, int i, int j, int k);
};

class SymmStencil : public BoundaryStencil<FlowField> {
 public:
  SymmStencil(const Parameters& parameters);

  void applyLeftWall(FlowField& flowField, int i, int j);
  void applyRightWall(FlowField& flowField, int i, int j);
  void applyBottomWall(FlowField& flowField, int i, int j);
  void applyTopWall(FlowField& flowField, int i, int j);

  void applyLeftWall(FlowField& flowField, int i, int j, int k);
  void applyRightWall(FlowField& flowField, int i, int j, int k);
  void applyBottomWall(FlowField& flowField, int i, int j, int k);
  void applyTopWall(FlowField& flowField, int i, int j, int k);
  void applyFrontWall(FlowField& flowField, int i, int j, int k);
  void applyBackWall(FlowField& flowField, int i, int j, int k);
};

class WallStencil : public BoundaryStencil<FlowField> {
 public:
  WallStencil(const Parameters& parameters);

  void applyLeftWall(FlowField& flowField, int i, int j);
  void applyRightWall(FlowField& flowField, int i, int j);
  void applyBottomWall(FlowField& flowField, int i, int j);
  void applyTopWall(FlowField& flowField, int i, int j);

  void applyLeftWall(FlowField& flowField, int i, int j, int k);
  void applyRightWall(FlowField& flowField, int i, int j, int k);
  void applyBottomWall(FlowField& flowField, int i, int j, int k);
  void applyTopWall(FlowField& flowField, int i, int j, int k);
  void applyFrontWall(FlowField& flowField, int i, int j, int k);
  void applyBackWall(FlowField& flowField, int i, int j, int k);
};
}
}
}

#endif
