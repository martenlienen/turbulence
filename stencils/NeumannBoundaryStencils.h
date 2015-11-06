#ifndef _NEUMANN_BOUNDARY_STENCIL_H_
#define _NEUMANN_BOUNDARY_STENCIL_H_

#include "../Parameters.h"
#include "../DataStructures.h"
#include "../FlowField.h"
#include "../Stencil.h"

class NeumannVelocityBoundaryStencil: public BoundaryStencil<FlowField> {

    public:

        /** Constructor
         * @param parameters Parameters of the simulation
         */
        NeumannVelocityBoundaryStencil(const Parameters & parameters);

        //@brief Functions for the 2D problem. Coordinates entered in alphabetical order.
        //@{
        void applyLeftWall   ( FlowField & flowField, int i, int j );
        void applyRightWall  ( FlowField & flowField, int i, int j );
        void applyBottomWall ( FlowField & flowField, int i, int j );
        void applyTopWall    ( FlowField & flowField, int i, int j );
        //@}

        //@brief Functions for the 3D problem. Coordinates entered in alphabetical order.
        //@{
        void applyLeftWall   ( FlowField & flowField, int i, int j, int k );
        void applyRightWall  ( FlowField & flowField, int i, int j, int k );
        void applyBottomWall ( FlowField & flowField, int i, int j, int k );
        void applyTopWall    ( FlowField & flowField, int i, int j, int k );
        void applyFrontWall  ( FlowField & flowField, int i, int j, int k );
        void applyBackWall   ( FlowField & flowField, int i, int j, int k );
        //@}
};

class NeumannFGHBoundaryStencil : public BoundaryStencil<FlowField> {

    public:

        /** Constructor
         * @param parameters Parameters of the simulation
         */
        NeumannFGHBoundaryStencil(const Parameters & parameters);

        //@brief Functions for the 2D problem. Coordinates entered in alphabetical order.
        //@{
        void applyLeftWall   ( FlowField & flowField, int i, int j );
        void applyRightWall  ( FlowField & flowField, int i, int j );
        void applyBottomWall ( FlowField & flowField, int i, int j );
        void applyTopWall    ( FlowField & flowField, int i, int j );
        //@}

        //@brief Functions for the 3D problem. Coordinates entered in alphabetical order.
        //@{
        void applyLeftWall   ( FlowField & flowField, int i, int j, int k );
        void applyRightWall  ( FlowField & flowField, int i, int j, int k );
        void applyBottomWall ( FlowField & flowField, int i, int j, int k );
        void applyTopWall    ( FlowField & flowField, int i, int j, int k );
        void applyFrontWall  ( FlowField & flowField, int i, int j, int k );
        void applyBackWall   ( FlowField & flowField, int i, int j, int k );
        //@}

};


#endif
