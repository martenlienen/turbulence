#ifndef _NSEOF_FLOWMODELS_LAMINAR_FLOW_FIELD_H_
#define _NSEOF_FLOWMODELS_LAMINAR_FLOW_FIELD_H_

#include "../../DataStructures.h"
#include "../../Parameters.h"
#include "../../FlowField.h"

namespace nseof {

namespace flowmodels {

namespace laminar {

/** Flow field
 *
 * Class intended to contain the state of the domain.
 */
class FlowField : public nseof::FlowField {
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
  FlowField(int Nx, int Ny, int Nz);
  FlowField(int Nx, int Ny);
  FlowField(Parameters const& parameters);

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
}
}

#endif
