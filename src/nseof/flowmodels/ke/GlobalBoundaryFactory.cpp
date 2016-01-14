#include "BFInputStencils.h"
#include "GlobalBoundaryFactory.h"

namespace nseof {

namespace flowmodels {

namespace ke {

GlobalBoundaryFactory::GlobalBoundaryFactory(Parameters& parameters)
    : nseof::GlobalBoundaryFactory(parameters) {
  if (parameters.simulation.scenario == "cavity") {
    // left
    _keStencils[0] = new WallStencil(parameters);

    // right
    _keStencils[1] = new WallStencil(parameters);
  } else {
    // left
    _keStencils[0] = new BFInputStencil(parameters);

    // right
    _keStencils[1] = new BFInputStencil(parameters);
  }
  // bottom
  _keStencils[2] = new WallStencil(parameters);

  // top
  if (parameters.simulation.scenario == "channel-symm") {
    _keStencils[3] = new SymmStencil(parameters);
  } else {
    _keStencils[3] = new WallStencil(parameters);
  }

  // front
  _keStencils[4] = new WallStencil(parameters);

  // back
  if (parameters.simulation.scenario == "channel-symm") {
    _keStencils[5] = new SymmStencil(parameters);
  } else {
    _keStencils[5] = new WallStencil(parameters);
  }
}

GlobalBoundaryFactory::~GlobalBoundaryFactory() {
  delete _keStencils[0];
  delete _keStencils[1];
  delete _keStencils[2];
  delete _keStencils[3];
  delete _keStencils[4];
  delete _keStencils[5];
}

GlobalBoundaryIterator<FlowField>
GlobalBoundaryFactory::getGlobalBoundaryIterator(FlowField& flowField) {
  if (_parameters.geometry.dim == 2) {
    return GlobalBoundaryIterator<FlowField>(
        flowField, _parameters, *(_keStencils[0]), *(_keStencils[1]),
        *(_keStencils[2]), *(_keStencils[3]), 1, 0);
  }
  return GlobalBoundaryIterator<FlowField>(
      flowField, _parameters, *(_keStencils[0]), *(_keStencils[1]),
      *(_keStencils[2]), *(_keStencils[3]), *(_keStencils[4]),
      *(_keStencils[5]), 1, 0);
}
}
}
}
