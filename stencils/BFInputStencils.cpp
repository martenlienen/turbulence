#include "BFInputStencils.h"

FLOAT computeVelocity3D (FlowField & flowField, int i, int j, int k, FLOAT stepSize,
                         const Parameters & parameters){
    const FLOAT posY = parameters.meshsize->getPosY(i,j,k);
    const FLOAT posZ = parameters.meshsize->getPosZ(i,j,k);
    const FLOAT dy   = parameters.meshsize->getDy(i,j,k);
    const FLOAT dz   = parameters.meshsize->getDz(i,j,k);

    if (posY+0.5*dy >= stepSize) {
        // Get the size of the inlet in Y. A 3 is subtracted because of the boundary cells
        const FLOAT inletYSize = parameters.geometry.lengthY - stepSize;
        const FLOAT inletZSize = parameters.geometry.lengthZ;

        const FLOAT y = posY+0.5*dy - stepSize;
        const FLOAT z = posZ+0.5*dz;

        return 36.0 * parameters.walls.vectorLeft[0] /
                      (inletZSize * inletZSize * inletYSize * inletYSize) *
                      y * (y - inletYSize) * z * (z - inletZSize);
    } else {
        return 0.0;
    }
}

FLOAT computeVelocity2D (FlowField & flowField, int i, int j, FLOAT stepSize,
                         const Parameters & parameters){
    const FLOAT posY = parameters.meshsize->getPosY(i,j);
    const FLOAT dy   = parameters.meshsize->getDy(i,j);

    if (posY+0.5*dy >= stepSize) {
        // Get the size of the inlet in Y. A 3 is subtracted because of the boundary cells
        const FLOAT inletYSize = parameters.geometry.lengthY - stepSize;

        const FLOAT y = posY+0.5*dy - stepSize;

        // DMITRIIS VERSION: for turbulence, please use: return parameters.walls.vectorLeft[0];
        return 6.0 * parameters.walls.vectorLeft[0] /
                     (inletYSize * inletYSize) * y * (inletYSize - y);
    } else {
        return 0.0;
    }
}

BFInputVelocityStencil::BFInputVelocityStencil (const Parameters & parameters) :
    BoundaryStencil<FlowField> (parameters),
    // Here, the obstacle size is set to zero if it was set as negative at the configuration
    _stepSize (parameters.bfStep.yRatio > 0.0 ? parameters.bfStep.yRatio*parameters.geometry.lengthY : 0.0)
{
       if (_parameters.geometry.dim==2){
               FLOAT posY = _parameters.meshsize->getPosY(0,0);
               FLOAT dy   = _parameters.meshsize->getDy(0,0);
               FLOAT nextDy = _parameters.meshsize->getDy(0,1);

               for (int j = 0; j < _parameters.geometry.sizeY-1; ++j) {
                       posY   = _parameters.meshsize->getPosY(0,j);
                       dy     = _parameters.meshsize->getDy(0,j);
                       nextDy = _parameters.meshsize->getDy(0,j+1);

                       // check if _stepSize is in this cell
                       if (posY + 0.5 * dy < _stepSize && _stepSize <=posY + dy + 0.5 * nextDy){
                               _stepSize = posY + dy;
                               break;
                       }
               }
       }
       else if(_parameters.geometry.dim==3){
               FLOAT posY = _parameters.meshsize->getPosY(0,0,0);
               FLOAT dy   = _parameters.meshsize->getDy(0,0,0);
               FLOAT nextDy = _parameters.meshsize->getDy(0,1,0);

               for (int j = 0; j < _parameters.geometry.sizeY-1; ++j) {
                       posY   = _parameters.meshsize->getPosY(0,j,0);
                       dy     = _parameters.meshsize->getDy(0,j,0);
                       nextDy = _parameters.meshsize->getDy(0,j+1,0);

                       // check if _stepSize is in this cell
                       if (posY + 0.5 * dy < _stepSize && _stepSize <=posY + dy + 0.5 * nextDy){
                               _stepSize = posY + dy;
                               break;
                       }
               }
       }
}

// Most of the functions are empty, and they shouldn't be called, assuming that the input is always
// located at the left.

void BFInputVelocityStencil::applyLeftWall   ( FlowField & flowField, int i, int j ){
    flowField.getVelocity().getVector(i,j)[0] =
        computeVelocity2D(flowField, i, j, _stepSize, _parameters);
    flowField.getVelocity().getVector(i,j)[1] = -flowField.getVelocity().getVector(i+1,j)[1];
}

void BFInputVelocityStencil::applyRightWall  ( FlowField & flowField, int i, int j ){}
void BFInputVelocityStencil::applyBottomWall ( FlowField & flowField, int i, int j ){}
void BFInputVelocityStencil::applyTopWall    ( FlowField & flowField, int i, int j ){}

void BFInputVelocityStencil::applyLeftWall   ( FlowField & flowField, int i, int j, int k ){
    flowField.getVelocity().getVector(i,j,k)[0] =
        computeVelocity3D(flowField, i, j, k, _stepSize, _parameters);
    flowField.getVelocity().getVector(i,j,k)[1] = -flowField.getVelocity().getVector(i+1,j,k)[1];
    flowField.getVelocity().getVector(i,j,k)[2] = -flowField.getVelocity().getVector(i+1,j,k)[2];
}

void BFInputVelocityStencil::applyRightWall  ( FlowField & flowField, int i, int j, int k ){}
void BFInputVelocityStencil::applyBottomWall ( FlowField & flowField, int i, int j, int k ){}
void BFInputVelocityStencil::applyTopWall    ( FlowField & flowField, int i, int j, int k ){}
void BFInputVelocityStencil::applyFrontWall  ( FlowField & flowField, int i, int j, int k ){}
void BFInputVelocityStencil::applyBackWall   ( FlowField & flowField, int i, int j, int k ){}


BFInputFGHStencil::BFInputFGHStencil(const Parameters & parameters) :
    BoundaryStencil<FlowField> (parameters),
    _stepSize (parameters.bfStep.yRatio> 0.0 ? parameters.bfStep.yRatio*parameters.geometry.lengthY : 0.0)
{}

void BFInputFGHStencil::applyLeftWall   ( FlowField & flowField, int i, int j ){
    flowField.getFGH().getVector(i,j)[0] =
        computeVelocity2D(flowField, i, j, _stepSize, _parameters);
}

void BFInputFGHStencil::applyRightWall  ( FlowField & flowField, int i, int j ){}
void BFInputFGHStencil::applyBottomWall ( FlowField & flowField, int i, int j ){}
void BFInputFGHStencil::applyTopWall    ( FlowField & flowField, int i, int j ){}

void BFInputFGHStencil::applyLeftWall  ( FlowField & flowField, int i, int j, int k ){
    flowField.getFGH().getVector(i,j,k)[0] =
        computeVelocity3D (flowField, i, j, k, _stepSize, _parameters);
}

void BFInputFGHStencil::applyRightWall  ( FlowField & flowField, int i, int j, int k ){}
void BFInputFGHStencil::applyBottomWall ( FlowField & flowField, int i, int j, int k ){}
void BFInputFGHStencil::applyTopWall    ( FlowField & flowField, int i, int j, int k ){}
void BFInputFGHStencil::applyFrontWall  ( FlowField & flowField, int i, int j, int k ){}
void BFInputFGHStencil::applyBackWall   ( FlowField & flowField, int i, int j, int k ){}
