#include "GlobalBoundaryFactory.h"

GlobalBoundaryFactory::GlobalBoundaryFactory(Parameters & parameters):
    _parameters(parameters){
    // The parameters will be modified, and therefore are not declared as constants.

    // All stencils are created, disregarding wether they will be used or not. This is less
    // complicated and doesn't seem that costly

    _periodic[0] = new PeriodicBoundaryVelocityStencil(parameters);
    _periodic[1] = new PeriodicBoundaryFGHStencil(parameters);

    _moving[0] = new MovingWallVelocityStencil(parameters);
    _moving[1] = new MovingWallFGHStencil(parameters);

    _outflow[0] = new NeumannVelocityBoundaryStencil(parameters);
    _outflow[1] = new NeumannFGHBoundaryStencil(parameters);

    _channelInput[0] = new BFInputVelocityStencil(parameters);
    _channelInput[1] = new BFInputFGHStencil(parameters);

    // Then, assign them according to the scenario
    std::string scenario = parameters.simulation.scenario;

    if (scenario == "cavity"){
        // Here, all is about setting the velocity at the boundaries
        for (int i = 0; i < 6; i++){
            _velocityStencils[i] = _moving[0];
            _FGHStencils[i] = _moving[1];
        }
        parameters.walls.typeLeft   = DIRICHLET;
        parameters.walls.typeRight  = DIRICHLET;
        parameters.walls.typeBottom = DIRICHLET;
        parameters.walls.typeTop    = DIRICHLET;
        parameters.walls.typeFront  = DIRICHLET;
        parameters.walls.typeBack   = DIRICHLET;
    } else if (scenario == "channel"){
        // To the left, we have the input
        _velocityStencils[0] = _channelInput[0];
        _FGHStencils[0] = _channelInput[1];

        // To the right, there is an outflow boundary
        _velocityStencils[1] = _outflow[0];
        _FGHStencils[1] = _outflow[1];

        // The other walls are moving walls
        for (int i = 2; i < 6; i++){
            _velocityStencils[i] = _moving[0];
            _FGHStencils[i] = _moving[1];
        }

        parameters.walls.typeLeft   = DIRICHLET;
        parameters.walls.typeRight  = NEUMANN;
        parameters.walls.typeBottom = DIRICHLET;
        parameters.walls.typeTop    = DIRICHLET;
        parameters.walls.typeFront  = DIRICHLET;
        parameters.walls.typeBack   = DIRICHLET;
    } else if ( scenario == "pressure-channel") {
      // we have Dirichlet conditions for pressure on both sides,
      // hence outflow conditions for the velocities
      _velocityStencils[0]=_outflow[0];
      _FGHStencils[0] = _outflow[1];

      // To the right, there is an outflow boundary
      _velocityStencils[1] = _outflow[0];
      _FGHStencils[1] = _outflow[1];

      // The other walls are moving walls
      for (int i = 2; i < 6; i++){
          _velocityStencils[i] = _moving[0];
          _FGHStencils[i] = _moving[1];
      }

      parameters.walls.typeLeft   = NEUMANN;
      parameters.walls.typeRight  = NEUMANN;
      parameters.walls.typeBottom = DIRICHLET;
      parameters.walls.typeTop    = DIRICHLET;
      parameters.walls.typeFront  = DIRICHLET;
      parameters.walls.typeBack   = DIRICHLET;
    } else if ( (scenario == "periodic-box") || (scenario== "taylor-green") ){
      for (int i = 0; i < 6; i++){
        _velocityStencils[i] = _periodic[0];
        _FGHStencils[i] = _periodic[1];
      }
      parameters.walls.typeLeft   = PERIODIC;
      parameters.walls.typeRight  = PERIODIC;
      parameters.walls.typeBottom = PERIODIC;
      parameters.walls.typeTop    = PERIODIC;
      parameters.walls.typeFront  = PERIODIC;
      parameters.walls.typeBack   = PERIODIC;
    } else {
        handleError(1, "Scenario not recognized");
    }
}

GlobalBoundaryFactory::~GlobalBoundaryFactory(){
    delete _moving[0];
    delete _moving[1];

    delete _periodic[0];
    delete _periodic[1];

    delete _outflow[0];
    delete _outflow[1];

    delete _channelInput[0];
    delete _channelInput[1];
}

GlobalBoundaryIterator<FlowField> GlobalBoundaryFactory::
    getGlobalBoundaryVelocityIterator(FlowField & flowField){
    if (_parameters.geometry.dim == 2){
        return GlobalBoundaryIterator<FlowField>(flowField, _parameters,
                                      *(_velocityStencils[0]), *(_velocityStencils[1]),
                                      *(_velocityStencils[2]), *(_velocityStencils[3]),
                                      1, 0);
    }
    return GlobalBoundaryIterator<FlowField>(flowField, _parameters,
                                  *(_velocityStencils[0]), *(_velocityStencils[1]),
                                  *(_velocityStencils[2]), *(_velocityStencils[3]),
                                  *(_velocityStencils[4]), *(_velocityStencils[5]),
                                  1, 0);
}

GlobalBoundaryIterator<FlowField> GlobalBoundaryFactory::
    getGlobalBoundaryFGHIterator(FlowField & flowField){
    if (_parameters.geometry.dim == 2){
        return GlobalBoundaryIterator<FlowField>(flowField, _parameters,
                                      *(_FGHStencils[0]), *(_FGHStencils[1]),
                                      *(_FGHStencils[2]), *(_FGHStencils[3]),
                                      1, 0);
    }
    return GlobalBoundaryIterator<FlowField>(flowField, _parameters,
                                  *(_FGHStencils[0]), *(_FGHStencils[1]),
                                  *(_FGHStencils[2]), *(_FGHStencils[3]),
                                  *(_FGHStencils[4]), *(_FGHStencils[5]),
                                  1, 0);
}
