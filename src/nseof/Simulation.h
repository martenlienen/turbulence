#ifndef _SIMULATION_H_
#define _SIMULATION_H_

#include "Parameters.h"

namespace nseof {

class Simulation {
 public:
  Simulation(Parameters &parameters) : _parameters(parameters) {}
  virtual ~Simulation() {}

  /** initialises the flow field according to the scenario */
  virtual void initializeFlowField() = 0;

  virtual void solveTimestep() = 0;

  virtual void plotVTK(int rank, int timeStep) = 0;

  virtual void serialize() = 0;
  virtual void deserialize() = 0;
  virtual void init() = 0;

 protected:
  Parameters &_parameters;
};
}

#endif  // _SIMULATION_H_
