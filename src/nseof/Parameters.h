#ifndef _PARAMETERS_H_
#define _PARAMETERS_H_

#include "Definitions.h"
#include <petscsys.h>
#include <string>
#include "Meshsize.h"

namespace nseof {

//! Classes for the parts of the parameters
//@{
class TimestepParameters {
 public:
  FLOAT dt;   //! Timestep
  FLOAT dtu;  //! Timestep (user)
  FLOAT tau;  //! Security factor
  FLOAT time = 0.0;
  FLOAT timeVTK = 0.0;
  int timeSteps = 0;
};

class SimulationParameters {
 public:
  FLOAT finalTime;       //! Final time for the simulation
  std::string type;      //! type of the simulation (DNS vs. turbulence)
  std::string scenario;  //! If channel or cavity, for example
  std::string uniform;
  std::string nulimiter;
};

class EnvironmentalParameters {
 public:
  // Gravity components
  FLOAT gx;
  FLOAT gy;
  FLOAT gz;
};

class FlowParameters {
 public:
  FLOAT Re;    //! Reynolds number
  FLOAT visc;  //! Reynolds number
  std::string type;
};

class kEpsilonParameters {
 public:
  FLOAT ce1;
  FLOAT ce2;
  FLOAT cmu;
  FLOAT sigmaK;
  FLOAT sigmaE;
  int model;
  FLOAT start;
  int adaptnrs;
  FLOAT adapterr;
};

class SolverParameters {
 public:
  FLOAT gamma;        //! Donor cell balance coefficient
  int maxIterations;  //! Maximum number of iterations in the linear solver
};

class GeometricParameters {
 public:
  // Dimensions
  int dim;

  // Number of cells
  int sizeX;
  int sizeY;
  int sizeZ;

  // Cell sizing
  FLOAT lengthX;
  FLOAT lengthY;
  FLOAT lengthZ;

  // meshsize type
  int meshsizeType;
  // for meshstretching
  int stretchX;
  int stretchY;
  int stretchZ;

  std::string obstacle;
};

class WallParameters {
 public:
  // Scalar value definition. Used to define the pressure, for example
  FLOAT scalarLeft;
  FLOAT scalarRight;
  FLOAT scalarBottom;
  FLOAT scalarTop;
  FLOAT scalarFront;
  FLOAT scalarBack;

  // Vector values at the boundaries, to define, for example, the velocities
  FLOAT vectorLeft[3];
  FLOAT vectorRight[3];
  FLOAT vectorBottom[3];
  FLOAT vectorTop[3];
  FLOAT vectorFront[3];
  FLOAT vectorBack[3];

  // Define how will the boundary behave
  BoundaryType typeLeft;
  BoundaryType typeRight;
  BoundaryType typeTop;
  BoundaryType typeBottom;
  BoundaryType typeFront;
  BoundaryType typeBack;
};

class VTKParameters {
 public:
  bool enabled;
  FLOAT interval;  //! Time interval for file printing
  int lowoffset;
  int highoffset;
  std::string prefix;  //! Output filename
  FLOAT start = 0.0;
};

class ParallelParameters {
 public:
  int rank;  //! Rank of the current processor

  int numProcessors[3];  //! Array with the number of processors in each
                         //! direction

  //@brief Ranks of the neighbors
  //@{
  int leftNb;
  int rightNb;
  int bottomNb;
  int topNb;
  int frontNb;
  int backNb;
  //@}

  int indices[3];      //! 3D indices to locate the array
  int localSize[3];    //! Size for the local flow field
  int firstCorner[3];  //! Position of the first element. Used for plotting

  PetscInt
      *sizes[3];  //! Arrays with the sizes of the blocks in each direction.
};

class BFStepParameters {
 public:
  FLOAT xRatio;
  FLOAT yRatio;
};

struct TimingParameters {
  bool enabled;
  std::string prefix;
};

struct RestartParameters {
  std::string in;
  std::string out;
};

//@}

/** A class to store and pass around the parameters
 */
class Parameters {
 public:
  Parameters() : meshsize(NULL) {}
  ~Parameters() {
    if (meshsize != NULL) {
      delete meshsize;
      meshsize = NULL;
    }
  }

  SimulationParameters simulation;
  TimestepParameters timestep;
  EnvironmentalParameters environment;
  FlowParameters flow;
  kEpsilonParameters kEpsilon;
  SolverParameters solver;
  GeometricParameters geometry;
  WallParameters walls;
  VTKParameters vtk;
  ParallelParameters parallel;
  BFStepParameters bfStep;
  TimingParameters timing;
  RestartParameters restart;
  // TODO WS2: include parameters for turbulence
  Meshsize *meshsize;
};
}

#endif
