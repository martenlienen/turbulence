#ifndef _NSEOF_FLOWMODELS_KE_NUT_STENCIL_H_
#define _NSEOF_FLOWMODELS_KE_NUT_STENCIL_H_

#include "../../Definitions.h"
#include "../../Parameters.h"
#include "../../Stencil.h"

#include "FlowField.h"

namespace nseof {

namespace flowmodels {

namespace ke {

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

class NutStencilA : public FieldStencil<FlowField> {
  virtual ~NutStencilA(){};

 public:
  FLOAT cmu;
  NutStencilA(const Parameters& parameters, FLOAT cmu);

  void apply(FlowField& flowField, int i, int j);

  void apply(FlowField& flowField, int i, int j, int k);
};
}
}
}

#endif
