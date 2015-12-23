#ifndef _STENCIL_NUT_H_
#define _STENCIL_NUT_H_

#include "../FlowFieldTurbA.h"
#include "../Stencil.h"
#include "../Parameters.h"

class NutStencil : public FieldStencil<FlowFieldTurbA> {
 public:
  NutStencil(const Parameters& parameters);

  ~NutStencil();

  void apply(FlowFieldTurbA& flowField, int i, int j);

  void apply(FlowFieldTurbA& flowField, int i, int j, int k);

 protected:
  FieldStencil<FlowFieldTurbA>* _fs;
};

class NutStencilL : public FieldStencil<FlowFieldTurbA> {
 public:
  NutStencilL(const Parameters& parameters);

  void apply(FlowFieldTurbA& flowField, int i, int j);

  void apply(FlowFieldTurbA& flowField, int i, int j, int k);
};

class NutStencilA : public FieldStencil<FlowFieldTurbA> {
 private:
  FLOAT _localVelocity[27 * 3];
  FLOAT _localMeshsize[27 * 3];

 public:
  NutStencilA(const Parameters& parameters);

  void apply(FlowFieldTurbA& flowField, int i, int j);

  void apply(FlowFieldTurbA& flowField, int i, int j, int k);
};

#endif
