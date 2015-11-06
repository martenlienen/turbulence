#ifndef _ITERATORS_H_
#define _ITERATORS_H_

#include "Stencil.h"
#include "Parameters.h"


/** Iterator class
 *
 * Applies operations to a flow field
 */
template<class FlowField>
class Iterator {
    protected:

        FlowField & _flowField;     //! Reference to the flow field
        const Parameters & _parameters; //! Rerence to the parameters

    public:

        /** Constructor for the iterator
         *
         * Creates an iterator, given a flow field and a stencil instance
         *
         * @param flowField Flow field with the state of the flow
         * @param stencil Stencil defining an operation to be applied
         */
        Iterator ( FlowField & flowfield, const Parameters& parameters ): _flowField(flowfield), _parameters(parameters){}

        /** Perform the stencil operation on inner, non-ghost cells
         */
        virtual void iterate () = 0;
};

template<class FlowField>
class FieldIterator : public Iterator<FlowField> {

    private:

        FieldStencil<FlowField> & _stencil;         //! Reference to a stencil

        //@brief Define the iteration domain to include more or less layers
        // Added since the ability to select the iteration domain provides more flexibility
        //@{
        const int _lowOffset;
        const int _highOffset;
        //@}

    public:

        FieldIterator (FlowField & flowField, const Parameters& parameters, FieldStencil<FlowField> & stencil,
                       int lowOffset = 0, int highOffset = 0);

        /** Volume iteration over the field.
         *
         * Volume iteration. The stencil will be applied to all cells in the domain plus the upper
         * boundaries. Lower boundaries are not included.
         */
        void iterate ();
};


template<class FlowField>
class GlobalBoundaryIterator : public Iterator<FlowField> {

    private:

        const int _lowOffset;
        const int _highOffset;

        // This iterator has a reference to a stencil for each side, and will call its methods
        BoundaryStencil<FlowField> & _leftWallStencil;   //! Stencil used on the left wall
        BoundaryStencil<FlowField> & _rightWallStencil;  //! Stencil used on the right wall
        BoundaryStencil<FlowField> & _bottomWallStencil; //! Stencil used on the bottom wall
        BoundaryStencil<FlowField> & _topWallStencil;    //! Stencil used on the top wall
        BoundaryStencil<FlowField> & _frontWallStencil;  //! Stencil used on the front wall
        BoundaryStencil<FlowField> & _backWallStencil;   //! Stencil used on the back wall

    public:

        /** Constructor for a single stencil in all faces
         * @param flowField The flowfield information
         * @param stencil Stencil for all faces
         */
        GlobalBoundaryIterator (FlowField & flowField, const Parameters & parameters,
                                BoundaryStencil<FlowField> & stencil,
                                int lowOffset = 0, int highOffset = 0);

        /** Constructor with different stencils for each face. For the 2D case.
         * @param flowField Flow field information
         * @param <some>WallStencil Stencil used on <some> wall
         */
        GlobalBoundaryIterator ( FlowField & flowField, const Parameters & parameters,
                                 BoundaryStencil<FlowField> & leftWallStencil,
                                 BoundaryStencil<FlowField> & rightWallStencil,
                                 BoundaryStencil<FlowField> & bottomWallStencil,
                                 BoundaryStencil<FlowField> & topWallStencil,
                                 int lowOffset = 0, int highOffset = 0);

        /** Constructor with different stencils for each face. For the 3D case.
         * @param flowField Flow field information
         * @param <some>WallStencil Stencil used on <some> wall
         */
        GlobalBoundaryIterator ( FlowField & flowField, const Parameters & parameters,
                                 BoundaryStencil<FlowField> & leftWallStencil,
                                 BoundaryStencil<FlowField> & rightWallStencil,
                                 BoundaryStencil<FlowField> & bottomWallStencil,
                                 BoundaryStencil<FlowField> & topWallStencil,
                                 BoundaryStencil<FlowField> & frontWallStencil,
                                 BoundaryStencil<FlowField> & backWallStencil,
                                 int lowOffset = 0, int highOffset = 0);

        /** Surface iterator
         *
         * Iterates on the boundary cells. Only upper corners and edges are iterated.
         */
        void iterate ();
};

template <class FlowField>
class ParallelBoundaryIterator : public Iterator<FlowField> {

    private:
        BoundaryStencil<FlowField> & _stencil;

        const int _lowOffset;
        const int _highOffset;

    public:
        ParallelBoundaryIterator(FlowField & flowField,
                                 const Parameters & parameters,
                                 BoundaryStencil<FlowField> & stencil,
                                 int lowOffset = 0, int highOffset = 0);
        void iterate();
};
#include "Iterators.cpph"

#endif
