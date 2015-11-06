#include "PeriodicBoundaryStencils.h"

PeriodicBoundaryVelocityStencil::PeriodicBoundaryVelocityStencil(const Parameters & parameters):
    BoundaryStencil<FlowField>(parameters) {}

// 2D problem

void PeriodicBoundaryVelocityStencil::applyLeftWall(FlowField & flowField, int i, int j){
    flowField.getVelocity().getVector(0, j)[0] =
        flowField.getVelocity().getVector(flowField.getNx(), j)[0];
    flowField.getVelocity().getVector(1, j)[1] =
        flowField.getVelocity().getVector(flowField.getNx()+1, j)[1];
}

void PeriodicBoundaryVelocityStencil::applyRightWall(FlowField & flowField, int i, int j){
    flowField.getVelocity().getVector(flowField.getNx()+2, j)[0] =
        flowField.getVelocity().getVector(2, j)[0];
    flowField.getVelocity().getVector(flowField.getNx()+2, j)[1] =
        flowField.getVelocity().getVector(2, j)[1];
}

void PeriodicBoundaryVelocityStencil::applyBottomWall(FlowField & flowField, int i, int j){
    flowField.getVelocity().getVector(i, 0)[1] =
        flowField.getVelocity().getVector(i, flowField.getNy())[1];
    flowField.getVelocity().getVector(i, 1)[0] =
        flowField.getVelocity().getVector(i, flowField.getNy()+1)[0];
}

void PeriodicBoundaryVelocityStencil::applyTopWall(FlowField & flowField, int i, int j){
    flowField.getVelocity().getVector(i, flowField.getNy()+2)[1] =
        flowField.getVelocity().getVector(i, 2)[1];
    flowField.getVelocity().getVector(i, flowField.getNy()+2)[0] =
        flowField.getVelocity().getVector(i, 2)[0];
}

// 3D Problem

void PeriodicBoundaryVelocityStencil::applyLeftWall(FlowField & flowField, int i, int j, int k){
    flowField.getVelocity().getVector(0, j, k)[0] =
        flowField.getVelocity().getVector(flowField.getNx(), j, k)[0];
    flowField.getVelocity().getVector(1, j, k)[1] =
        flowField.getVelocity().getVector(flowField.getNx()+1, j, k)[1];
    flowField.getVelocity().getVector(1, j, k)[2] =
        flowField.getVelocity().getVector(flowField.getNx()+1, j, k)[2];
}

void PeriodicBoundaryVelocityStencil::applyRightWall(FlowField & flowField, int i, int j, int k){
    flowField.getVelocity().getVector(flowField.getNx()+2, j, k)[0] =
        flowField.getVelocity().getVector(2, j, k)[0];
    flowField.getVelocity().getVector(flowField.getNx()+2, j, k)[1] =
        flowField.getVelocity().getVector(2, j, k)[1];
    flowField.getVelocity().getVector(flowField.getNx()+2, j, k)[2] =
        flowField.getVelocity().getVector(2, j, k)[2];
}

void PeriodicBoundaryVelocityStencil::applyBottomWall(FlowField & flowField, int i, int j, int k){
    flowField.getVelocity().getVector(i, 1, k)[0] =
        flowField.getVelocity().getVector(i, flowField.getNy()+1, k)[0];
    flowField.getVelocity().getVector(i, 0, k)[1] =
        flowField.getVelocity().getVector(i, flowField.getNy(), k)[1];
    flowField.getVelocity().getVector(i, 1, k)[2] =
        flowField.getVelocity().getVector(i, flowField.getNy()+1, k)[2];
}

void PeriodicBoundaryVelocityStencil::applyTopWall(FlowField & flowField, int i, int j, int k){
    flowField.getVelocity().getVector(i, flowField.getNy()+2, k)[0] =
        flowField.getVelocity().getVector(i, 2, k)[0];
    flowField.getVelocity().getVector(i, flowField.getNy()+2, k)[1] =
        flowField.getVelocity().getVector(i, 2, k)[1];
    flowField.getVelocity().getVector(i, flowField.getNy()+2, k)[2] =
        flowField.getVelocity().getVector(i, 2, k)[2];
}
void PeriodicBoundaryVelocityStencil::applyFrontWall(FlowField & flowField, int i, int j, int k){
    flowField.getVelocity().getVector(i, j, 1)[0] =
        flowField.getVelocity().getVector(i, j, flowField.getNz()+1)[0];
    flowField.getVelocity().getVector(i, j, 1)[1] =
        flowField.getVelocity().getVector(i, j, flowField.getNz()+1)[1];
    flowField.getVelocity().getVector(i, j, 0)[2] =
        flowField.getVelocity().getVector(i, j, flowField.getNz())[2];
}

void PeriodicBoundaryVelocityStencil::applyBackWall(FlowField & flowField, int i, int j, int k){
    flowField.getVelocity().getVector(i, j, flowField.getNz()+2)[0] =
        flowField.getVelocity().getVector(i, j, 2)[0];
    flowField.getVelocity().getVector(i, j, flowField.getNz()+2)[1] =
        flowField.getVelocity().getVector(i, j, 2)[1];
    flowField.getVelocity().getVector(i, j, flowField.getNz()+2)[2] =
        flowField.getVelocity().getVector(i, j, 2)[2];
}


PeriodicBoundaryFGHStencil::PeriodicBoundaryFGHStencil(const Parameters & parameters):
    BoundaryStencil<FlowField>(parameters) {}

// 2D problem
void PeriodicBoundaryFGHStencil::applyLeftWall(FlowField & flowField, int i, int j){}
void PeriodicBoundaryFGHStencil::applyRightWall(FlowField & flowField, int i, int j){}
void PeriodicBoundaryFGHStencil::applyBottomWall(FlowField & flowField, int i, int j){}
void PeriodicBoundaryFGHStencil::applyTopWall(FlowField & flowField, int i, int j){}

// 3D Problem
void PeriodicBoundaryFGHStencil::applyLeftWall(FlowField & flowField, int i, int j, int k){}
void PeriodicBoundaryFGHStencil::applyRightWall(FlowField & flowField, int i, int j, int k){}
void PeriodicBoundaryFGHStencil::applyBottomWall(FlowField & flowField, int i, int j, int k){}
void PeriodicBoundaryFGHStencil::applyTopWall(FlowField & flowField, int i, int j, int k){}
void PeriodicBoundaryFGHStencil::applyFrontWall(FlowField & flowField, int i, int j, int k){}
void PeriodicBoundaryFGHStencil::applyBackWall(FlowField & flowField, int i, int j, int k){}
