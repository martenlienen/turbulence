#ifndef _NSEOF_FLOWMODELS_ALGEBRAIC_NUTSTENCIL_H_
#define _NSEOF_FLOWMODELS_ALGEBRAIC_NUTSTENCIL_H_

#include "../../Stencil.h"
#include "../../Parameters.h"

#include "FlowField.h"

namespace nseof {

namespace flowmodels {

namespace algebraic {

class NutStencil : public FieldStencil<FlowField> {
 public:
  NutStencil(const Parameters& parameters);

  ~NutStencil();

  void apply(FlowField& flowField, int i, int j);

  void apply(FlowField& flowField, int i, int j, int k);

 protected:
  FieldStencil<FlowField>* _fs;
};

class NutStencilL : public FieldStencil<FlowField> {
 public:
  NutStencilL(const Parameters& parameters);

  void apply(FlowField& flowField, int i, int j);

  void apply(FlowField& flowField, int i, int j, int k);
};

class limiter {
 public:
  limiter() {}
  virtual ~limiter() {}
  virtual void limit(const Parameters& parameters, FLOAT& lm, int i, int j) = 0;

  virtual void limit(const Parameters& parameters, FLOAT& lm, int i, int j,
                     int k) = 0;
};

class NutStencilA : public FieldStencil<FlowField> {
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

  void apply(FlowField& flowField, int i, int j);

  void apply(FlowField& flowField, int i, int j, int k);
};
}
}
}

#endif
