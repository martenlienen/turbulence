#ifndef _STENCIL_NUT_H_
#define _STENCIL_NUT_H_

#include "../FlowFieldTurbA.h"
#include "../Stencil.h"
#include "../Parameters.h"

namespace nseof {

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

class limiter {
 public:
  limiter() {}
  virtual ~limiter() {}
  virtual void limit(const Parameters& parameters, FLOAT& lm, int i, int j) = 0;

  virtual void limit(const Parameters& parameters, FLOAT& lm, int i, int j,
                     int k) = 0;
};

class NutStencilA : public FieldStencil<FlowFieldTurbA> {
 private:
  FLOAT _localVelocity[27 * 3];
  FLOAT _localMeshsize[27 * 3];
  limiter* l;

  virtual ~NutStencilA() { delete l; };

  void computeNUT2D(int i, int j, const FLOAT* const localVelocity,
                    const FLOAT* const localMeshsize,
                    const Parameters& parameters, const FLOAT& h, FLOAT& nu,
                    FLOAT& flu, FLOAT& lmm);

  void computeNUT3D(int i, int j, int k, const FLOAT* const localVelocity,
                    const FLOAT* const localMeshsize,
                    const Parameters& parameters, const FLOAT& h, FLOAT& nu,
                    FLOAT& flu, FLOAT& lmm);

 public:
  NutStencilA(const Parameters& parameters);

  void apply(FlowFieldTurbA& flowField, int i, int j);

  void apply(FlowFieldTurbA& flowField, int i, int j, int k);
};
}

#endif
