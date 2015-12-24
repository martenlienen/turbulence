#ifndef _FLOW_FIELD_TURB_A_H_
#define _FLOW_FIELD_TURB_A_H_

#include "DataStructures.h"
#include "Parameters.h"
#include "FlowField.h"
#include "FlowFieldLaminar.h"

namespace nseof {

/** Flow field for algebraic turbulent model:
 *    - it is a flow field
 *    - has a laminar flow field
 *    - has additional scalar fields
 *        (vortex viscosity, wall distance, velocity fluctuation, mixing length)
 *
 */
class FlowFieldTurbA : public FlowField {
 private:
  FlowFieldLaminar _flowField;
  ScalarField _nu;
  ScalarField _h;
  ScalarField _u;
  ScalarField _lm;

 public:
  FlowFieldTurbA(const Parameters& parameters);

  int getNx() const;

  int getNy() const;

  int getNz() const;

  int getCellsX() const;

  int getCellsY() const;

  int getCellsZ() const;

  ScalarField& getPressure();

  VectorField& getVelocity();

  IntScalarField& getFlags();

  VectorField& getFGH();

  ScalarField& getRHS();

  void getPressureAndVelocity(FLOAT& pressure, FLOAT* const velocity, int i,
                              int j);
  void getPressureAndVelocity(FLOAT& pressure, FLOAT* const velocity, int i,
                              int j, int k);

  // vortex viscosity
  FLOAT& getNu(int i, int j);
  FLOAT& getNu(int i, int j, int k);

  // wall distance
  FLOAT& getH(int i, int j);
  FLOAT& getH(int i, int j, int k);

  // velocity fluctuation
  FLOAT& getU(int i, int j);
  FLOAT& getU(int i, int j, int k);

  // mixing length
  FLOAT& getLm(int i, int j);
  FLOAT& getLm(int i, int j, int k);
};
}

#endif
