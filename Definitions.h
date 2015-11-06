#ifndef _DEFINITIONS_H_
#define _DEFINITIONS_H_

#include <stdlib.h>
#include <stdio.h>
#include <iostream>
#include <float.h>
#include <petscksp.h>

// Datatype for the type of data stored in the structures
#ifdef USE_SINGLE_PRECISION
#define MY_MPI_FLOAT MPI_FLOAT
#define MY_FLOAT_MAX FLT_MAX
typedef float FLOAT;
#else
#define MY_MPI_FLOAT MPI_DOUBLE
#define MY_FLOAT_MAX DBL_MAX
typedef double FLOAT;
#endif


// Types of boundary
enum BoundaryType {
    DIRICHLET,        // Also used for the case of non-moving wall
    PERIODIC,
    PARALLEL_BOUNDARY,
    NEUMANN
};

const int OBSTACLE_SELF =   1<<0;
const int OBSTACLE_LEFT =   1<<1;
const int OBSTACLE_RIGHT =  1<<2;
const int OBSTACLE_BOTTOM = 1<<3;
const int OBSTACLE_TOP =    1<<4;
const int OBSTACLE_FRONT =  1<<5;
const int OBSTACLE_BACK =   1<<6;


// An assertion sending back a message
#ifdef DEBUG
#define assertion(boolean)\
    if (!(boolean)){std::cerr << "Assertion failed: " << #boolean\
        << " In file " << __FILE__ << ", line " << __LINE__ << "."\
        << std::endl; exit(2);}
#else
#define assertion(boolean)
#endif

#endif


// An error handler, at least to check memory allocation. Active also outside
// debug

#define handleError(code,message)\
    { std::cerr << "ERROR in file " << __FILE__ << ", line " << __LINE__\
              << std::endl << #message << std::endl;\
    exit ( code ); }
