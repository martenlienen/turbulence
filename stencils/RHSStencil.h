#ifndef _RHS_STENCIL_H_
#define _RHS_STENCIL_H_

#include "../Stencil.h"
#include "../Parameters.h"
#include "../FlowField.h"

/** Field stencil to compute the right hand side of the pressure equation.
 */
class RHSStencil : public FieldStencil<FlowField> {

    public:

        /** Constructs and instance of the RHS Stencil
         *
         * @param parameters Parameters of the flow
         */
        RHSStencil ( const Parameters & parameters );

        /** Apply the stencil in 2D
         * @param flowField Flow field to work on
         * @param i Position in the X direction
         * @param j Position in the Y direction
         */
        void apply ( FlowField & flowField, int i, int j );

        /** Apply the stencil in 3D
         * @param flowField Flow field to work on
         * @param i Position in the X direction
         * @param j Position in the Y direction
         * @param k Position in the Z direction
         */
        void apply ( FlowField & flowField, int i, int j, int k );
};

#endif
