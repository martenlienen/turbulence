#ifndef _DATA_STRUCTURES_H_
#define _DATA_STRUCTURES_H_

#include "Definitions.h"

/** Storage of a scalar field
 *
 * Parent of storage classes. Contains the data pointer and sizes in each
 * dimension.
 */
template<class DataType>
class Field
{
    protected:

        //! Pointer to the data array
        DataType* _data;

        const int _size_x; //! Size of the field in x direction, including ghost layers
        const int _size_y; //! Size of the field in y direction, including ghost layers
        const int _size_z; //! Size of the field in z direction, including ghost layers
        const int _components; //! Number of components per position
        const int _size;       //! Total size of the data array


    public:

        /** Constructor for the field
         *
         * General constructor. Takes the three arguments even if the matrix is
         * two dimensional. Doesn't allocate memory, only sets the dimensions
         *
         * @param Nx Number of cells in the x direction
         * @param Ny Number of cells in the y direction
         * @param Nz Number of cells in the z direction
         */
        Field ( int Nx, int Ny, int Nz, int components ):
          _size_x ( Nx ), _size_y ( Ny ), _size_z ( Nz ), _components ( components ),
          _size ( components * Nx * Ny * Nz ) {

          _data = new DataType [ _size ];

          // Check if the data was allocated successfully
          if (_data == 0){
              handleError ( 1, "Unable to allocate memory" );
          }
        }

        /** Destructor
         *
         * Deallocates the data for all the fields
         */
        virtual ~Field (){
          if (_data != NULL){
              delete[] _data;
              _data = NULL;
          }
        }

        /** Returns the number of cells in the x direction
         *
         * @return The size in the x direction
         */
        int getNx () const { return _size_x;}

        /** Returns the number of cells in the y direction
         *
         * @return The size in the y direction
         */
        int getNy () const { return _size_y;}

        /** Returns the number of cells in the z direction
         *
         * @return The size in the z direction
         */
        int getNz () const { return _size_z;}

        /** Index to array position mapper
         *
         * Index mapper. Converts the given index to the corresponding position
         * in the array.
         *
         * @param i x index
         * @param j y index
         * @param k z index
         *
         * @return Position in the array
         */
        int index2array ( int i, int j, int k = 0 ) const {

          // Check that the indexes are within bounds
          assertion ( ( i < _size_x ) && ( j < _size_y ) && ( k < _size_z ) );
          assertion ( ( i >= 0 ) && ( j >= 0 ) && ( k >= 0 )  );

          return  _components * ( i + ( j * _size_x ) + ( k * _size_x * _size_y ) );
        }

};


/** Scalar field representation
 *
 * Stores a scalar field of floats. Derived from Field.
 */
class ScalarField: public Field<FLOAT>
{
    private:
        /** initializes the data of the scalar field by zeros. */
        void initialize();

    public:

        /** 2D scalar field constructor.
         *
         * Sets the size of the data array and allocates data for the 2D field
         *
         * @param Nx Number of cells in direction x
         * @param Ny Number of cells in direction y
         */
        ScalarField ( int Nx, int Ny );

        /** 3D scalar field constructor.
         *
         * Sets the size of the data array and allocates data for the 3D field
         *
         * @param Nx Number of cells in direction x
         * @param Ny Number of cells in direction y
         * @param Nz Number of cells in direction z
         */
        ScalarField ( int Nx, int Ny, int Nz );

        /** Acces to element in scalar field
         *
         * Returns a reference to an element of the scalar field, so that it
         * can be read and written.
         * @param i x index
         * @param j y index
         * @param k z index. Not required for arrays of dimension two.
         */
        FLOAT & getScalar ( int i, int j, int k = 0 );

        /** Prints the contents of the field
         *
         * Shows the content of the scalar field by printing them to stdout. Only clear if the
         * matrix is small enough. Used for debugging.
         *
         * @param title A label for the printed matrix
         */
        void show (const std::string title = "");

};


/** Vector field representation
 *
 * Stores a vector field of floats. Derived from Field.
 */
class VectorField: public Field<FLOAT>
{
    private:
         /** initializes the data of the scalar field by zeros. */
         void initialize();

    public:

        /** 2D Vector field constructor.
         *
         * Sets the size of the data array and allocates data for the 2D field
         *
         * @param Nx Number of cells in direction x
         * @param Ny Number of cells in direction y
         * @param Nz Number of cells in direction z
         */
        VectorField ( int Nx, int Ny );

        /** 3D Vector field constructor.
         *
         * Sets the size of the data array and allocates data for the 3D field
         *
         * @param Nx Number of cells in direction x
         * @param Ny Number of cells in direction y
         * @param Nz Number of cells in direction z
         */
        VectorField ( int Nx, int Ny, int Nz );

        /** Non constant acces to an element in the vector field
         *
         * Returns a pointer to the position in the array that can be used to
         * modify it.
         *
         * @param i x index
         * @param j y index
         * @param k z index
         */
        FLOAT* getVector ( int i, int j, int k = 0 );

        /** Prints the contents of the field
         *
         * Shows the content of the first two components of a vector field by printing them to
         * stdout. Only clear if the matrix is small enough. Used for debugging.
         *
         * @param title A label for the printed matrix
         */
        void show (const std::string title = "");
};


/** Integer field
 * 
 * Integer field with one value per position. Intended to represent flag
 * fields. Implemented because templates are undesirable at this point
 */
class IntScalarField : public Field<int> {

    private:

        /** Initialize the array with zeros
         *
         * Had to be redeclared here because the object of the parent class
         * writes floats and caused a memory access violation.
         */
        void initialize ();

    public:

        /** 2D constructor
         *
         * @param Nx Size in the x direction
         * @param Ny Size in the Y direction
         */
        IntScalarField ( int Nx, int Ny );

        /** 3D constructor
         *
         * @param Nx Size in the x direction
         * @param Ny Size in the Y direction
         * @param Nz SIze in the Z direction
         */
        IntScalarField ( int Nx, int Ny, int Nz );

        /** Access field values
         *
         * Returns a reference to the element with the given index
         *
         * @param i X index
         * @param j Y index
         * @param k Z index
         */
        int & getValue ( int i, int j, int k = 0 );

        void show (const std::string title = "");
};

#endif
