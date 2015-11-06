#include "PetscParallelConfiguration.h"

PetscParallelConfiguration::PetscParallelConfiguration(Parameters & parameters):
    _parameters(parameters) {

    // Obtain the rank of the current processor
    int rank, nproc;
    MPI_Comm_rank(PETSC_COMM_WORLD, &rank);
    MPI_Comm_size(PETSC_COMM_WORLD, &nproc);

    parameters.parallel.rank = rank;

    // Obtain the position of this subdomain, and locate its neighbors.
    createIndices();
    locateNeighbors();
    computeSizes();

    int nprocFromFile = _parameters.parallel.numProcessors[0] *
                        _parameters.parallel.numProcessors[1];

    if (parameters.geometry.dim == 3) {
        nprocFromFile *= _parameters.parallel.numProcessors[2];
    }

    if (nproc != nprocFromFile){
        handleError(1, "The number of processors specified in the configuration file doesn't match the communicator");
    }
}


PetscParallelConfiguration::~PetscParallelConfiguration(){
    freeSizes();
}


void PetscParallelConfiguration::locateNeighbors(){

    int i = _parameters.parallel.indices[0];
    int j = _parameters.parallel.indices[1];
    int k = _parameters.parallel.indices[2];

    if (_parameters.geometry.dim == 2){
        _parameters.parallel.leftNb   = computeRankFromIndices(i-1, j, 0);
        _parameters.parallel.rightNb  = computeRankFromIndices(i+1, j, 0);
        _parameters.parallel.bottomNb = computeRankFromIndices(i, j-1, 0);
        _parameters.parallel.topNb    = computeRankFromIndices(i, j+1, 0);

        // The following two are not used in this case
        _parameters.parallel.frontNb = MPI_PROC_NULL;
        _parameters.parallel.backNb  = MPI_PROC_NULL;
    } else {
        _parameters.parallel.leftNb   = computeRankFromIndices(i-1, j, k);
        _parameters.parallel.rightNb  = computeRankFromIndices(i+1, j, k);
        _parameters.parallel.bottomNb = computeRankFromIndices(i, j-1, k);
        _parameters.parallel.topNb    = computeRankFromIndices(i, j+1, k);
        _parameters.parallel.frontNb  = computeRankFromIndices(i, j, k-1);
        _parameters.parallel.backNb   = computeRankFromIndices(i, j, k+1);
    }

    // If periodic boundaries declared, let the process itself deal with it, without communication
    int rank;
    MPI_Comm_rank(PETSC_COMM_WORLD, &rank);
}


void PetscParallelConfiguration::createIndices(){
    int & rank = _parameters.parallel.rank;
    _parameters.parallel.indices[0] = rank % _parameters.parallel.numProcessors[0];
    _parameters.parallel.indices[1] = (rank / _parameters.parallel.numProcessors[0]) %
                                      _parameters.parallel.numProcessors[1];
    _parameters.parallel.indices[2] = rank / (_parameters.parallel.numProcessors[0] *
                                              _parameters.parallel.numProcessors[1]);
}


int PetscParallelConfiguration::computeRankFromIndices(int i, int j, int k) const {

    if( i < 0 || i >= _parameters.parallel.numProcessors[0] ||
        j < 0 || j >= _parameters.parallel.numProcessors[1] ||
        k < 0 || k >= _parameters.parallel.numProcessors[2] ){
        return MPI_PROC_NULL;
    }
    int nrank = i + j*_parameters.parallel.numProcessors[0];
    if (_parameters.geometry.dim == 3){
        nrank += k * _parameters.parallel.numProcessors[0] * _parameters.parallel.numProcessors[1];
    }
    return nrank;
}


void PetscParallelConfiguration::computeSizes(){

    int dim = _parameters.geometry.dim;

    for (int i = 0; i < dim; i++){
        _parameters.parallel.sizes[i] = new PetscInt[_parameters.parallel.numProcessors[i]];
    }

    int geometrySizes[3];
    geometrySizes[0] = _parameters.geometry.sizeX;
    geometrySizes[1] = _parameters.geometry.sizeY;
    geometrySizes[2] = _parameters.geometry.sizeZ;

    for (int i = 0; i < dim; i++){
        for (int j = 0; j < _parameters.parallel.numProcessors[i]; j++){
            _parameters.parallel.sizes[i][j] =
                geometrySizes[i] / _parameters.parallel.numProcessors[i];
            if (j < geometrySizes[i] % _parameters.parallel.numProcessors[i]){
                _parameters.parallel.sizes[i][j] ++;
            }
        }
    }

    // Locate the position of the first element of the subdomain. Useful for plotting later on.
    for (int i = 0; i < dim; i++){
        _parameters.parallel.firstCorner[i] = 0;
        for (int j = 0; j < _parameters.parallel.indices[i]; j++){
            _parameters.parallel.firstCorner[i] += _parameters.parallel.sizes[i][j];
        }
    }
    if (dim == 2){
        _parameters.parallel.firstCorner[2] = 0;
    }

    // Select the local sizes from the already computed sizes
    for (int i = 0; i < dim; i++){
        _parameters.parallel.localSize[i] =
            _parameters.parallel.sizes[i][_parameters.parallel.indices[i]];
    }

    // If the domain lies on an edge, add one to that direction, for the artificial external
    // pressures in the PETSc solver
    for (int i = 0; i < dim; i++){
        _parameters.parallel.sizes[i][0] ++;
        _parameters.parallel.sizes[i][_parameters.parallel.numProcessors[i]-1] ++;
    }
}


void PetscParallelConfiguration::freeSizes(){

    int dim = _parameters.geometry.dim;

    for (int i = 0; i < dim; i++){
        delete[] _parameters.parallel.sizes[i];
    }
}
