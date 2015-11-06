#include "RHSStencil.h"

RHSStencil::RHSStencil ( const Parameters & parameters ) : FieldStencil<FlowField> ( parameters ) {}


void RHSStencil::apply ( FlowField & flowField, int i, int j ) {
    flowField.getRHS().getScalar (i, j) = 1.0 / _parameters.timestep.dt *
        ( ( flowField.getFGH().getVector(i, j)[0] - flowField.getFGH().getVector(i-1, j)[0] )
          / _parameters.meshsize->getDx(i,j)
        + ( flowField.getFGH().getVector(i, j)[1] - flowField.getFGH().getVector(i, j-1)[1] )
          / _parameters.meshsize->getDy(i,j) );
}


void RHSStencil::apply ( FlowField & flowField, int i, int j, int k ) {
    flowField.getRHS().getScalar (i, j, k) = 1.0 / _parameters.timestep.dt *
        ( (flowField.getFGH().getVector(i, j, k)[0] - flowField.getFGH().getVector(i-1, j, k)[0])
          / _parameters.meshsize->getDx(i,j,k)
        + (flowField.getFGH().getVector(i, j, k)[1] - flowField.getFGH().getVector(i, j-1, k)[1])
          / _parameters.meshsize->getDy(i,j,k)
        + (flowField.getFGH().getVector(i, j, k)[2] - flowField.getFGH().getVector(i, j, k-1)[2])
          / _parameters.meshsize->getDz(i,j,k) );
}
