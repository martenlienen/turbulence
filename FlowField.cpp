#include "FlowField.h"

FlowField::FlowField ( int Nx, int Ny ) :
    _size_x ( Nx ), _size_y ( Ny ), _size_z ( 1 ),
    _cellsX (Nx+3), _cellsY(Ny+3), _cellsZ(1),

    // Pressure field doesn't need to have an extra layer, but this allows to address the same
    // positions with the same iterator for both pressures and velocities.
    _pressure ( ScalarField ( Nx + 3, Ny + 3 ) ),
    _velocity ( VectorField ( Nx + 3, Ny + 3 ) ), _flags ( IntScalarField ( Nx + 3, Ny + 3 ) ),
    _FGH ( VectorField ( Nx + 3, Ny + 3 ) ), _RHS ( ScalarField (Nx + 3, Ny + 3) ) {

    assertion ( Nx > 0 );
    assertion ( Ny > 0 );
}


FlowField::FlowField ( int Nx, int Ny, int Nz ) :
    _size_x ( Nx ), _size_y ( Ny ), _size_z ( Nz ),
    _cellsX (Nx+3), _cellsY(Ny+3), _cellsZ(Nz+3),
    _pressure ( ScalarField ( Nx + 3, Ny + 3, Nz + 3 ) ),
    _velocity  ( VectorField ( Nx + 3, Ny + 3, Nz + 3 ) ),
    _flags  ( IntScalarField ( Nx + 3, Ny + 3, Nz +3 ) ),
    _FGH ( VectorField ( Nx + 3, Ny + 3, Nz + 3 ) ),
    _RHS ( ScalarField ( Nx + 3, Ny + 3, Nz + 3 ) ) {


    // Check that the provided data makes sense
    assertion ( Nx > 0 );
    assertion ( Ny > 0 );
    assertion ( Nz > 0 );
}


FlowField::FlowField (const Parameters & parameters):
    _size_x(parameters.parallel.localSize[0]),
    _size_y(parameters.parallel.localSize[1]),
    _size_z(parameters.parallel.localSize[2]),
    _cellsX(_size_x+3),
    _cellsY(_size_y+3),
    _cellsZ(_size_z+3),
    // Probably far from the best way to write this
    _pressure(parameters.geometry.dim==2?ScalarField(_size_x + 3, _size_y + 3):
                      ScalarField(_size_x + 3, _size_y + 3, _size_z + 3)),
    _velocity(parameters.geometry.dim==2?VectorField(_size_x + 3, _size_y + 3):
                      VectorField(_size_x + 3, _size_y + 3, _size_z + 3)),
    _flags(parameters.geometry.dim==2?IntScalarField(_size_x + 3, _size_y + 3):
                   IntScalarField(_size_x + 3, _size_y + 3, _size_z + 3)),
    _FGH(parameters.geometry.dim==2?VectorField(_size_x + 3, _size_y + 3):
                 VectorField(_size_x + 3, _size_y + 3, _size_z + 3)),
    _RHS(parameters.geometry.dim==2?ScalarField(_size_x + 3, _size_y + 3):
                 ScalarField(_size_x + 3, _size_y + 3, _size_z + 3))
{ }


int FlowField::getNx () const{
    return _size_x;
}


int FlowField::getNy () const{
    return _size_y;
}


int FlowField::getNz () const{
    return _size_z;
}


int FlowField::getCellsX() const{
    return _cellsX;
}

int FlowField::getCellsY() const{
    return _cellsY;
}

int FlowField::getCellsZ() const{
    return _cellsZ;
}


ScalarField & FlowField::getPressure () {
    return _pressure;
}


VectorField & FlowField::getVelocity () {
    return _velocity;
}


IntScalarField & FlowField::getFlags () {
    return _flags;
}


VectorField & FlowField::getFGH () {
    return _FGH;
}


ScalarField & FlowField::getRHS () {
    return _RHS;
}

void FlowField::getPressureAndVelocity(FLOAT &pressure, FLOAT* const velocity,  int i, int j){
    FLOAT * v_here = getVelocity().getVector(i, j);
    FLOAT * v_left = getVelocity().getVector(i-1, j);
    FLOAT * v_down = getVelocity().getVector(i, j-1);

    velocity[0] = ( v_here[0] + v_left[0] ) / 2;
    velocity[1] = ( v_here[1] + v_down[1] ) / 2;

    pressure = getPressure().getScalar(i,j);
}

void FlowField::getPressureAndVelocity(FLOAT &pressure, FLOAT* const velocity, int i, int j, int k){
    FLOAT * v_here = getVelocity().getVector(i, j, k);
    FLOAT * v_left = getVelocity().getVector(i-1, j, k);
    FLOAT * v_down = getVelocity().getVector(i, j-1, k);
    FLOAT * v_back = getVelocity().getVector(i, j, k-1);

    velocity[0] = ( v_here[0] + v_left[0] ) / 2;
    velocity[1] = ( v_here[1] + v_down[1] ) / 2;
    velocity[2] = ( v_here[2] + v_back[2] ) / 2;

    pressure = getPressure().getScalar(i,j,k);
}
