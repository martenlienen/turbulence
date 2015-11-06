#ifndef _VELOCITY_STENCIL_H_
#define _VELOCITY_STENCIL_H_

#include "../Stencil.h"
#include "../Parameters.h"
#include "../FlowField.h"

/** Stencil to compute the velocity once the pressure has been found.
 */
class VelocityStencil : public FieldStencil<FlowField> {

    public:

        /** Constructor
         * @param parameters Parameters of the problem
         */
        VelocityStencil(const Parameters & parameters);

        /** Apply the stencil in 2D
         * @param flowField Flow field information
         * @param i Position in the X direction
         * @param j Position in the Y direction
         */
        void apply ( FlowField & flowField, int i, int j );

        /** Apply the stencil in 3D
         * @param flowField Flow field information
         * @param i Position in the X direction
         * @param j Position in the Y direction
         * @param k Position in the Z direction
         */
        void apply ( FlowField & flowField, int i, int j, int k );
};

#endif
