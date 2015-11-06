#ifndef _PERIODIC_BOUNDARY_STENCIL_H_
#define _PERIODIC_BOUNDARY_STENCIL_H_

#include "../Stencil.h"
#include "../Parameters.h"
#include "../FlowField.h"

/** Stencil to set periodic boundary conditions for velocity
 */
class PeriodicBoundaryVelocityStencil: public BoundaryStencil<FlowField> {

    public:

        /** Constructor
         * @param parameters Parameters of the simulation
         */
        PeriodicBoundaryVelocityStencil(const Parameters & parameters);

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


/** Stencil to set periodic boundary conditions for velocity for FGH. Since there are no operations
 * in FGH, this stencil does nothing.
 */
class PeriodicBoundaryFGHStencil: public BoundaryStencil<FlowField> {

    public:

        /** Constructor
         * @param parameters Parameters of the simulation
         */
        PeriodicBoundaryFGHStencil(const Parameters & parameters);

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
