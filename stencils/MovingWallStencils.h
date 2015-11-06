#ifndef _MOVING_WALL_STENCIL_H_
#define _MOVING_WALL_STENCIL_H_

#include "../Stencil.h"
#include "../Parameters.h"
#include "../FlowField.h"

/** Boundary stencil to set constant velocities at the faces.
 *
 * The values for the velocities are stored in the parameters.
 */
class MovingWallVelocityStencil: public BoundaryStencil<FlowField> {

    public:

        /** Constructor
         *
         * @param parameters Parameters of the problem
         */
        MovingWallVelocityStencil ( const Parameters & parameters );

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

/** Boundary stencil to set correct values of FGH at the faces.
 *
 * The values for the velocities are stored in the parameters.
 */
class MovingWallFGHStencil: public BoundaryStencil<FlowField> {

    public:

        /** Constructor
         *
         * @param parameters Parameters of the problem
         */
        MovingWallFGHStencil ( const Parameters & parameters );

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
