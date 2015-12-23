#ifndef _FLOW_FIELD_H_
#define _FLOW_FIELD_H_

#include "DataStructures.h"
#include "Parameters.h"

/** Flow field
 *
 * Class intended to contain the state of the domain.
 */
class FlowField {
 public:
  //  /** Constructor for the 2D flow field
  //   *
  //   * Constructor for the flow field. Allocates all the fields and sets
  //   * the sizes. Currently, this contructor is only used for testing
  //   purposes.
  //   *
  //   * @param Nx Size of the fuild domain (non-ghost cells), in the X
  //   direction
  //   * @param Ny Size of the fuild domain (non-ghost cells), in the Y
  //   direction
  //   */
  //  FlowField(int Nx, int Ny, int Nz);
  //
  //  /** Constructor for the 3D flow field
  //   *
  //   * Constructor for the flow field. Allocates all the fields and sets
  //   * the sizes. Currently, this contructor is only used for testing
  //   purposes.
  //   *
  //   * @param Nx Size of the fuild domain (non-ghost cells), in the X
  //   direction
  //   * @param Ny Size of the fuild domain (non-ghost cells), in the Y
  //   direction
  //   * @param Nz Size of the fuild domain (non-ghost cells), in the Z
  //   direction
  //   */
  //  FlowField(int Nx, int Ny);
  //
  //  /** Constructs a field from parameters object
  //   *
  //   * Constructs a field from a parameters object, so that it dimensionality
  //   can
  //   * be defined in
  //   * the configuration file.
  //   *
  //   * @param parameters Parameters object with geometric information
  //   */
  //  FlowField(const Parameters& parameters);

  virtual ~FlowField() {}

  /** Obtain size in the X direction
   *
   * @return Number of cells in the X direction
   */
  virtual int getNx() const = 0;

  /** Obtain size in the Y direction
   *
   * @return Number of cells in the Y direction
   */
  virtual int getNy() const = 0;

  /** Obtain size in the Z direction
   *
   * @return Number of cells in the Z direction
   */
  virtual int getNz() const = 0;

  virtual int getCellsX() const = 0;
  virtual int getCellsY() const = 0;
  virtual int getCellsZ() const = 0;

  /** Get pressure field
   * @return Reference to pressure field
   */
  virtual ScalarField& getPressure() = 0;

  /** Get velocity field
   * @return Reference to velocity field
   */
  virtual VectorField& getVelocity() = 0;

  /** Get flag field
   * @return Reference to flag field
   */
  virtual IntScalarField& getFlags() = 0;

  /** Get the field with the F, G, and H  abbreviations
   * @return Multi-component field with F, G and H
   */
  virtual VectorField& getFGH() = 0;

  /** Get the right hand side for the pressure linear solver
   * @return Scalar field with the right hand side
   */
  virtual ScalarField& getRHS() = 0;

  virtual void getPressureAndVelocity(FLOAT& pressure, FLOAT* const velocity,
                                      int i, int j) = 0;
  virtual void getPressureAndVelocity(FLOAT& pressure, FLOAT* const velocity,
                                      int i, int j, int k) = 0;
};

#endif
