#include "NeumannBoundaryStencils.h"

NeumannVelocityBoundaryStencil::NeumannVelocityBoundaryStencil(const Parameters & parameters):
    BoundaryStencil<FlowField>(parameters) {}


void NeumannVelocityBoundaryStencil::applyLeftWall   ( FlowField & flowField, int i, int j ){
    flowField.getVelocity().getVector(i-1,j)[0] = flowField.getVelocity().getVector(i,j)[0];
    flowField.getVelocity().getVector(i,j)[1] = flowField.getVelocity().getVector(i+1,j)[1];
}

void NeumannVelocityBoundaryStencil::applyRightWall  ( FlowField & flowField, int i, int j ){
    flowField.getVelocity().getVector(i,j)[0] = flowField.getVelocity().getVector(i-1,j)[0];
    flowField.getVelocity().getVector(i,j)[1] = flowField.getVelocity().getVector(i-1,j)[1];
}

void NeumannVelocityBoundaryStencil::applyBottomWall ( FlowField & flowField, int i, int j ){
    flowField.getVelocity().getVector(i,j)[0] = flowField.getVelocity().getVector(i,j+1)[0];
    flowField.getVelocity().getVector(i,j-1)[1] = flowField.getVelocity().getVector(i,j)[1];
}

void NeumannVelocityBoundaryStencil::applyTopWall    ( FlowField & flowField, int i, int j ){
    flowField.getVelocity().getVector(i,j)[0] = flowField.getVelocity().getVector(i,j-1)[0];
    flowField.getVelocity().getVector(i,j)[1] = flowField.getVelocity().getVector(i,j-1)[1];
}

void NeumannVelocityBoundaryStencil::applyLeftWall   ( FlowField & flowField, int i, int j, int k ){
    flowField.getVelocity().getVector(i-1,j,k)[0] = flowField.getVelocity().getVector(i,j,k)[0];
    flowField.getVelocity().getVector(i,j,k)[1] = flowField.getVelocity().getVector(i+1,j,k)[1];
    flowField.getVelocity().getVector(i,j,k)[2] = flowField.getVelocity().getVector(i+1,j,k)[2];
}

void NeumannVelocityBoundaryStencil::applyRightWall  ( FlowField & flowField, int i, int j, int k ){
    flowField.getVelocity().getVector(i,j,k)[0] = flowField.getVelocity().getVector(i-1,j,k)[0];
    flowField.getVelocity().getVector(i,j,k)[1] = flowField.getVelocity().getVector(i-1,j,k)[1];
    flowField.getVelocity().getVector(i,j,k)[2] = flowField.getVelocity().getVector(i-1,j,k)[2];
}

void NeumannVelocityBoundaryStencil::applyBottomWall ( FlowField & flowField, int i, int j, int k ){
    flowField.getVelocity().getVector(i,j,k)[0] = flowField.getVelocity().getVector(i,j+1,k)[0];
    flowField.getVelocity().getVector(i,j-1,k)[1] = flowField.getVelocity().getVector(i,j,k)[1];
    flowField.getVelocity().getVector(i,j,k)[2] = flowField.getVelocity().getVector(i,j+1,k)[2];
}

void NeumannVelocityBoundaryStencil::applyTopWall    ( FlowField & flowField, int i, int j, int k ){
    flowField.getVelocity().getVector(i,j,k)[0] = flowField.getVelocity().getVector(i,j-1,k)[0];
    flowField.getVelocity().getVector(i,j,k)[1] = flowField.getVelocity().getVector(i,j-1,k)[1];
    flowField.getVelocity().getVector(i,j,k)[2] = flowField.getVelocity().getVector(i,j-1,k)[2];
}

void NeumannVelocityBoundaryStencil::applyFrontWall  ( FlowField & flowField, int i, int j, int k ){
    flowField.getVelocity().getVector(i,j,k)[0] = flowField.getVelocity().getVector(i,j,k+1)[0];
    flowField.getVelocity().getVector(i,j,k)[1] = flowField.getVelocity().getVector(i,j,k+1)[1];
    flowField.getVelocity().getVector(i,j,k-1)[2] = flowField.getVelocity().getVector(i,j,k)[2];
}

void NeumannVelocityBoundaryStencil::applyBackWall   ( FlowField & flowField, int i, int j, int k ){
    flowField.getVelocity().getVector(i,j,k)[0] = flowField.getVelocity().getVector(i,j,k-1)[0];
    flowField.getVelocity().getVector(i,j,k)[1] = flowField.getVelocity().getVector(i,j,k-1)[1];
    flowField.getVelocity().getVector(i,j,k)[2] = flowField.getVelocity().getVector(i,j,k-1)[2];
}


NeumannFGHBoundaryStencil::NeumannFGHBoundaryStencil (const Parameters & parameters) :
    BoundaryStencil<FlowField> (parameters)
{}

// These are left empty. The right values should be computed by the FGH body stencil

void NeumannFGHBoundaryStencil::applyLeftWall ( FlowField & flowField, int i, int j ){}
void NeumannFGHBoundaryStencil::applyRightWall ( FlowField & flowField, int i, int j ){}
void NeumannFGHBoundaryStencil::applyBottomWall ( FlowField & flowField, int i, int j ){}
void NeumannFGHBoundaryStencil::applyTopWall ( FlowField & flowField, int i, int j ){}


// 3D stencils

void NeumannFGHBoundaryStencil::applyLeftWall ( FlowField & flowField, int i, int j, int k ){}
void NeumannFGHBoundaryStencil::applyRightWall ( FlowField & flowField, int i, int j , int k ){}
void NeumannFGHBoundaryStencil::applyBottomWall ( FlowField & flowField, int i, int j, int k ){}
void NeumannFGHBoundaryStencil::applyTopWall ( FlowField & flowField, int i, int j, int k ){}
void NeumannFGHBoundaryStencil::applyFrontWall ( FlowField & flowField, int i, int j, int k ){}
void NeumannFGHBoundaryStencil::applyBackWall ( FlowField & flowField, int i, int j, int k ){}
