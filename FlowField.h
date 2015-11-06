#ifndef _FLOW_FIELD_H_
#define _FLOW_FIELD_H_

#include "DataStructures.h"
#include "Parameters.h"

/** Flow field
 *
 * Class intended to contain the state of the domain.
 */
class FlowField {

    private:

        const int _size_x; //! Size in the X direction
        const int _size_y; //! Size in the Y direction
        const int _size_z; //! Size in the Z direction

        const int _cellsX;
        const int _cellsY;
        const int _cellsZ;

        ScalarField _pressure; //! Scalar field representing the pressure
        VectorField _velocity; //! Multicomponent field representing velocity
        IntScalarField _flags; //! Integer field for the flags

        VectorField _FGH;

        ScalarField _RHS;      //! Right hand side for the Poisson equation

    public:

        /** Constructor for the 2D flow field
         *
         * Constructor for the flow field. Allocates all the fields and sets
         * the sizes. Currently, this contructor is only used for testing purposes.
         *
         * @param Nx Size of the fuild domain (non-ghost cells), in the X direction
         * @param Ny Size of the fuild domain (non-ghost cells), in the Y direction
         */
        FlowField (int Nx, int Ny, int Nz);

        /** Constructor for the 3D flow field
         *
         * Constructor for the flow field. Allocates all the fields and sets
         * the sizes. Currently, this contructor is only used for testing purposes.
         *
         * @param Nx Size of the fuild domain (non-ghost cells), in the X direction
         * @param Ny Size of the fuild domain (non-ghost cells), in the Y direction
         * @param Nz Size of the fuild domain (non-ghost cells), in the Z direction
         */
        FlowField ( int Nx, int Ny );

        /** Constructs a field from parameters object
         *
         * Constructs a field from a parameters object, so that it dimensionality can be defined in
         * the configuration file.
         *
         * @param parameters Parameters object with geometric information
         */
        FlowField (const Parameters & parameters);

        /** Obtain size in the X direction
         *
         * @return Number of cells in the X direction
         */
        int getNx () const;

        /** Obtain size in the Y direction
         *
         * @return Number of cells in the Y direction
         */
        int getNy () const;

        /** Obtain size in the Z direction
         *
         * @return Number of cells in the Z direction
         */
        int getNz () const;

        int getCellsX() const;
        int getCellsY() const;
        int getCellsZ() const;

        /** Get pressure field
         * @return Reference to pressure field
         */
        ScalarField & getPressure ();

        /** Get velocity field
         * @return Reference to velocity field
         */
        VectorField & getVelocity ();

        /** Get flag field
         * @return Reference to flag field
         */
        IntScalarField & getFlags ();

        /** Get the field with the F, G, and H  abbreviations
         * @return Multi-component field with F, G and H
         */
        VectorField & getFGH ();

        /** Get the right hand side for the pressure linear solver
         * @return Scalar field with the right hand side
         */
        ScalarField & getRHS ();

        void getPressureAndVelocity(FLOAT &pressure, FLOAT* const velocity, int i, int j);
        void getPressureAndVelocity(FLOAT &pressure, FLOAT* const velocity, int i, int j, int k);
};

#endif
