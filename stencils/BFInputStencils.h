#ifndef _BF_INPUT_VELOCITY_STENCIL_H_
#define _BF_INPUT_VELOCITY_STENCIL_H_

#include "../Stencil.h"
#include "../Parameters.h"
#include "../FlowField.h"

/**
 * A stencil to set the input velocity in channel flows. Only implements the applyLeftWall(...) methods.
 */
class BFInputVelocityStencil : public BoundaryStencil<FlowField> {

    public:
        BFInputVelocityStencil (const Parameters & parameters);

        void applyLeftWall   ( FlowField & flowField, int i, int j );
        void applyRightWall  ( FlowField & flowField, int i, int j );
        void applyBottomWall ( FlowField & flowField, int i, int j );
        void applyTopWall    ( FlowField & flowField, int i, int j );

        void applyLeftWall   ( FlowField & flowField, int i, int j, int k );
        void applyRightWall  ( FlowField & flowField, int i, int j, int k );
        void applyBottomWall ( FlowField & flowField, int i, int j, int k );
        void applyTopWall    ( FlowField & flowField, int i, int j, int k );
        void applyFrontWall  ( FlowField & flowField, int i, int j, int k );
        void applyBackWall   ( FlowField & flowField, int i, int j, int k );

    private:
        FLOAT _stepSize; //! fixes the size of the step. If zero, is just channel flow
};

/** FGH stencil which corresponds to one-sided Dirichlet conditions (i.e. the channel flow profile).
 *  Only implements the applyLeftWall(...) methods.
 */
class BFInputFGHStencil : public BoundaryStencil<FlowField> {

    public:
        BFInputFGHStencil (const Parameters & parameters);

        void applyLeftWall   ( FlowField & flowField, int i, int j );
        void applyRightWall  ( FlowField & flowField, int i, int j );
        void applyBottomWall ( FlowField & flowField, int i, int j );
        void applyTopWall    ( FlowField & flowField, int i, int j );

        void applyLeftWall   ( FlowField & flowField, int i, int j, int k );
        void applyRightWall  ( FlowField & flowField, int i, int j, int k );
        void applyBottomWall ( FlowField & flowField, int i, int j, int k );
        void applyTopWall    ( FlowField & flowField, int i, int j, int k );
        void applyFrontWall  ( FlowField & flowField, int i, int j, int k );
        void applyBackWall   ( FlowField & flowField, int i, int j, int k );

    private:
        FLOAT _stepSize; //! fixes the size of the step. If zero, is just channel flow
};

#endif
