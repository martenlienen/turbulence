#include "BFStepInitStencil.h"

BFStepInitStencil::BFStepInitStencil (const Parameters & parameters) :
    FieldStencil<FlowField> (parameters),
    xLimit(parameters.bfStep.xRatio*parameters.geometry.lengthX),
    yLimit(parameters.bfStep.yRatio*parameters.geometry.lengthY)
{}

void BFStepInitStencil::apply(FlowField & flowField, int i, int j){
    IntScalarField & flags = flowField.getFlags();
    const FLOAT posX = _parameters.meshsize->getPosX(i,j);
    const FLOAT posY = _parameters.meshsize->getPosY(i,j);
    const FLOAT dx   = _parameters.meshsize->getDx(i,j);
    const FLOAT dy   = _parameters.meshsize->getDy(i,j);
    const FLOAT nextDx   = _parameters.meshsize->getDx(i+1,j);
    const FLOAT nextDy   = _parameters.meshsize->getDy(i,j+1);
    const FLOAT lastDx   = _parameters.meshsize->getDx(i-1,j);
    const FLOAT lastDy   = _parameters.meshsize->getDy(i,j-1);

    if (posX+0.5*dx < xLimit && posY+0.5*dy < yLimit){
       flags.getValue(i, j) = OBSTACLE_SELF;
    }
    if (posX-0.5*lastDx < xLimit && posY+0.5*dy < yLimit){
        flags.getValue(i, j) += OBSTACLE_LEFT;
    }
    if (posX+dx+0.5*nextDx < xLimit && posY+0.5*dy < yLimit){
        flags.getValue(i, j) += OBSTACLE_RIGHT;
    }
    if (posX+0.5*dx < xLimit && posY-0.5*lastDy < yLimit){
        flags.getValue(i, j) += OBSTACLE_BOTTOM;
    }
    if (posX+0.5*dx < xLimit && posY+dy+0.5*nextDy < yLimit){
        flags.getValue(i, j) += OBSTACLE_TOP;
    }
}

void BFStepInitStencil::apply(FlowField & flowField, int i, int j, int k){
    IntScalarField & flags = flowField.getFlags();
    const FLOAT posX = _parameters.meshsize->getPosX(i,j,k);
    const FLOAT posY = _parameters.meshsize->getPosY(i,j,k);
    const FLOAT dx   = _parameters.meshsize->getDx(i,j,k);
    const FLOAT dy   = _parameters.meshsize->getDy(i,j,k);
    const FLOAT nextDx   = _parameters.meshsize->getDx(i+1,j,k);
    const FLOAT nextDy   = _parameters.meshsize->getDy(i,j+1,k);
    const FLOAT lastDx   = _parameters.meshsize->getDx(i-1,j,k);
    const FLOAT lastDy   = _parameters.meshsize->getDy(i,j-1,k);

    if (posX+0.5*dx < xLimit && posY+0.5*dy < yLimit){
        flags.getValue(i, j, k) = OBSTACLE_SELF;

        // The obstacle is 2D, so we can say the following
        flags.getValue(i, j, k) += OBSTACLE_FRONT;
        flags.getValue(i, j, k) += OBSTACLE_BACK;
    }
    if (posX-0.5*lastDx < xLimit && posY+0.5*dy < yLimit){
        flags.getValue(i, j, k) += OBSTACLE_LEFT;
    }
    if (posX+dx+0.5*nextDx < xLimit && posY+0.5*dy < yLimit){
        flags.getValue(i, j, k) += OBSTACLE_RIGHT;
    }
    if (posX+0.5*dx < xLimit && posY-0.5*lastDy < yLimit){
        flags.getValue(i, j, k) += OBSTACLE_BOTTOM;
    }
    if (posX+0.5*dx < xLimit && posY+dx+0.5*nextDy < yLimit){
        flags.getValue(i, j, k) += OBSTACLE_TOP;
    }
}
