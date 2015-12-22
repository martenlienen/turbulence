#ifndef _FLOW_FIELD_TURB_A_H_
#define _FLOW_FIELD_TURB_A_H_

#include "DataStructures.h"
#include "Parameters.h"
#include "FlowField.h"
#include "FlowFieldLaminar.h"

/** Flow field
 *
 * Class intended to contain the state of the domain.
 */
class FlowFieldTurbA : public FlowField {
 private:
  FlowFieldLaminar _flowField;
  ScalarField _nu;
  ScalarField _h;

 public:
  /** Constructs a field from parameters object
   *
   * Constructs a field from a parameters object, so that it dimensionality can
   * be defined in
   * the configuration file.
   *
   * @param parameters Parameters object with geometric information
   */
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

  FLOAT& getNu(int i, int j);
  FLOAT& getNu(int i, int j, int k);

  FLOAT& getH(int i, int j);
  FLOAT& getH(int i, int j, int k);
};

#endif
