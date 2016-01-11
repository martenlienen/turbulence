#ifndef _NSEOF_FLOWMODEL_KE_FLOW_FIELD_H_
#define _NSEOF_FLOWMODEL_KE_FLOW_FIELD_H_

#include "../../DataStructures.h"
#include "../../Parameters.h"

#include "../turbulent/FlowField.h"

namespace nseof {

namespace flowmodels {

namespace ke {

/**
 * Flow field for a k-epsilon turbulence model
 */
class FlowField : public nseof::flowmodels::turbulent::FlowField {
 private:
  ScalarField _fmu;
  ScalarField _f1;
  ScalarField _f2;
  ScalarField _f3;
  ScalarField _D;
  ScalarField _E;
  ScalarField _epsilon;
  ScalarField _tke;
  ScalarField _rhs_epsilon;
  ScalarField _rhs_tke;
  ScalarField _sijsij;

 public:
  FlowField(const Parameters& parameters);

  ScalarField& getFmu();
  ScalarField& getF1();
  ScalarField& getF2();
  ScalarField& getF3();
  ScalarField& getD();
  ScalarField& getE();
  ScalarField& getEpsilon();
  ScalarField& getTke();
  ScalarField& getRHSEpsilon();
  ScalarField& getRHSTke();
  ScalarField& getsijsij();

  // mixing length
  FLOAT& getFmu(int i, int j);
  FLOAT& getFmu(int i, int j, int k);
  FLOAT& getF1(int i, int j);
  FLOAT& getF1(int i, int j, int k);
  FLOAT& getF2(int i, int j);
  FLOAT& getF2(int i, int j, int k);
  FLOAT& getF3(int i, int j);
  FLOAT& getF3(int i, int j, int k);
  FLOAT& getD(int i, int j);
  FLOAT& getD(int i, int j, int k);
  FLOAT& getE(int i, int j);
  FLOAT& getE(int i, int j, int k);
  FLOAT& getEpsilon(int i, int j);
  FLOAT& getEpsilon(int i, int j, int k);
  FLOAT& getTke(int i, int j);
  FLOAT& getTke(int i, int j, int k);
  FLOAT& getRHSEpsilon(int i, int j);
  FLOAT& getRHSEpsilon(int i, int j, int k);
  FLOAT& getRHSTke(int i, int j);
  FLOAT& getRHSTke(int i, int j, int k);
  FLOAT& getsijsij(int i, int j);
  FLOAT& getsijsij(int i, int j, int k);
};
}
}
}

#endif
