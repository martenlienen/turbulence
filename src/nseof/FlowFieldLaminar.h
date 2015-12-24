#ifndef _FLOW_FIELD_LAMINAR_H_
#define _FLOW_FIELD_LAMINAR_H_

#include "DataStructures.h"
#include "Parameters.h"
#include "FlowField.h"

namespace nseof {

/** Flow field
 *
 * Class intended to contain the state of the domain.
 */
class FlowFieldLaminar : public FlowField {
 private:
  const int _size_x;  //! Size in the X direction
  const int _size_y;  //! Size in the Y direction
  const int _size_z;  //! Size in the Z direction

  const int _cellsX;
  const int _cellsY;
  const int _cellsZ;

  ScalarField _pressure;  //! Scalar field representing the pressure
  VectorField _velocity;  //! Multicomponent field representing velocity
  IntScalarField _flags;  //! Integer field for the flags

  VectorField _FGH;

  ScalarField _RHS;  //! Right hand side for the Poisson equation

 public:
  /** Constructs a field from parameters object
   *
   * Constructs a field from a parameters object, so that it dimensionality can
   * be defined in
   * the configuration file.
   *
   * @param parameters Parameters object with geometric information
   */
  FlowFieldLaminar(int Nx, int Ny, int Nz);
  FlowFieldLaminar(int Nx, int Ny);
  FlowFieldLaminar(const Parameters& parameters);

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
};
}

#endif
