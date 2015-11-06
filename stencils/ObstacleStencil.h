#ifndef _OBSTACLE_STENCIL_H_
#define _OBSTACLE_STENCIL_H_

#include "../Stencil.h"
#include "../Parameters.h"
#include "../FlowField.h"

/** Compute all velocities on obstacle cells. This has been taken out of
 * Velocity stencil to circumvent any race conditions
 **/
class ObstacleStencil : public FieldStencil<FlowField> {
	public:
		/* Constructor
		 * @param parameters Parameters of the problem
		 */
		ObstacleStencil(const Parameters & parameters);


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
