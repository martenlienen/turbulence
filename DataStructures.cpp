#include "DataStructures.h"


// Functions for scalar and vector fields -------------------------------------


// Functions specific to the scalar field -------------------------------------


ScalarField::ScalarField ( int Nx, int Ny ): Field<FLOAT>( Nx, Ny, 1, 1 ) {
    initialize();
}


ScalarField::ScalarField ( int Nx, int Ny, int Nz ): Field<FLOAT> ( Nx, Ny, Nz, 1 ) {
    initialize();
}


FLOAT& ScalarField::getScalar ( int i, int j, int k ) {
    return _data [ index2array ( i, j, k ) ];
}

void ScalarField::show(const std::string title){
    std::cout << std::endl << "--- " << title << " ---" << std::endl;
    for (int k = 0; k < _size_z; k++){
        for (int j = _size_y-1; j > -1; j--){
            for (int i = 0; i < _size_x; i++){
                std::cout << getScalar(i,j,k) << "\t";
            }
            std::cout << std::endl;
        }
        std::cout << std::endl;
    }
}

void ScalarField::initialize () {
    for ( int i = 0; i < _size; i++ ){
        _data[i] = 0.0;
    }
}



// Functions related to the vector field --------------------------------------

VectorField::VectorField ( int Nx, int Ny ): Field<FLOAT> ( Nx, Ny, 1,  2 ) {
    initialize();
}


VectorField::VectorField ( int Nx, int Ny, int Nz ): Field<FLOAT> ( Nx, Ny, Nz, 3 ) {
    initialize();
}


FLOAT* VectorField::getVector ( int i, int j, int k ) {
    return &_data[index2array( i, j, k )];
}


void VectorField::show(const std::string title){
    std::cout << std::endl << "--- " << title << " ---" << std::endl;
    std::cout << "Component 1" << std::endl;
    for (int k = 0; k < _size_z; k++){
        for (int j = _size_y-1; j > -1; j--){
            for (int i = 0; i < _size_x; i++){
                std::cout << getVector(i,j,k)[0] << "\t";
            }
            std::cout << std::endl;
        }
        std::cout << std::endl;
    }

    std::cout << "Component 2" << std::endl;
    for (int k = 0; k < _size_z; k++){
        for (int j = _size_y-1; j > -1; j--){
            for (int i = 0; i < _size_x; i++){
                std::cout << getVector(i,j,k)[1] << "\t";
            }
            std::cout << std::endl;
        }
        std::cout << std::endl;
    }
}

void VectorField::initialize () {
    for ( int i = 0; i < _size; i++ ){
        _data[i] = 0.0;
    }
}

//----------------------------------------------------------------------------
// Stuff related to the integer scalar field

IntScalarField::IntScalarField ( int Nx, int Ny ) : Field<int> ( Nx, Ny, 1, 1 ) {
  initialize();
}


IntScalarField::IntScalarField ( int Nx, int Ny, int Nz ) : Field<int> ( Nx, Ny, Nz, 1 ) {
  initialize();
}


void IntScalarField::initialize () {
    for ( int i = 0; i < _size; i++ ){
        _data[i] = 0;
    }
}


int & IntScalarField::getValue ( int i, int j, int k ) {
    return _data[ index2array (i, j, k) ];
}


void IntScalarField::show(const std::string title){
    std::cout << std::endl << "--- " << title << " ---" << std::endl;
    for (int k = 0; k < _size_z; k++){
        for (int j = _size_y-1; j > -1; j--){
            for (int i = 0; i < _size_x; i++){
                std::cout << getValue(i,j,k) << "\t";
            }
            std::cout << std::endl;
        }
        std::cout << std::endl;
    }
}
