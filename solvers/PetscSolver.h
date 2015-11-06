#ifndef _PETSC_SOLVER_H_
#define _PETSC_SOLVER_H_

#include <petscksp.h>
#include <petscdm.h>
#include <petscdmda.h>
#include "../FlowField.h"
#include "../DataStructures.h"
#include "../Parameters.h"
#include "../LinearSolver.h"

const unsigned char LEFT_WALL_BIT   = 1<<0;
const unsigned char RIGHT_WALL_BIT  = 1<<1;
const unsigned char BOTTOM_WALL_BIT = 1<<2;
const unsigned char TOP_WALL_BIT    = 1<<3;
const unsigned char FRONT_WALL_BIT  = 1<<4;
const unsigned char BACK_WALL_BIT   = 1<<5;


/** A class to encapsulate information the Petsc builder functions
 *  Petsc used so called context objects to give information to its routines.
 *  We need them to pass the flow field and the parameters in a single argument.
 */
class PetscUserCtx{
    private:
        Parameters & _parameters;  //! Reference to parameters
        FlowField & _flowField;  //! Reference to the flow field

        int *_limitsX, *_limitsY, *_limitsZ;

        int _rank;

    public:

        /** Constructor
         *
         * @param parameters A parameters instance
         * @param flowField The corresponding flow field
         */
        PetscUserCtx(Parameters & parameters, FlowField & flowField);

        /** Returns the parameters */
        Parameters & getParameters();

        /** Returns the flow field */
        FlowField & getFlowField();

        void setLimits(int *limitsX, int *limitsY, int *limitsZ);
        void getLimits(int ** limitsX, int ** limitsY, int ** limitsZ);

        void setRank (int rank);
        int getRank() const;

        unsigned char setAsBoundary;    // if set as boundary in the linear system. Use bits
        int displacement[6];            // Displacements for the boundary treatment

};


class PetscSolver : public LinearSolver{

    private:
        Vec _x;  //! Petsc vectors for solution and RHS
        DM _da;  //! Topology manager
        KSP _ksp;  //! Solver context
	PC _pc;  //! Preconditioner

        PetscUserCtx _ctx;        // Capsule for Petsc builders

        // Indices for filling the matrices and right hand side
        int _limitsX[2], _limitsY[2], _limitsZ[2];

        PetscInt _firstX, _lengthX, _firstY, _lengthY, _firstZ, _lengthZ;

        // Additional variables used to determine where to write back the results
        int _offsetX, _offsetY, _offsetZ;

    public:

        /** Constructor */
        PetscSolver(FlowField & flowField, Parameters & parameters);

        /** Uses petsc to solve the linear system for the pressure */
        void solve();

        /** Returns the grid */
        const DM & getGrid() const;

	/** Reinit the matrix so that it uses the right flag field */
	void reInitMatrix();

};

#endif
