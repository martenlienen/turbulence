#include "ObstacleStencil.h"

ObstacleStencil::ObstacleStencil(const Parameters & parameters) :
	FieldStencil<FlowField> (parameters) {}

void ObstacleStencil::apply ( FlowField & flowField, int i, int j){
	const int obstacle = flowField.getFlags().getValue(i, j);
	VectorField & velocity = flowField.getVelocity();    

	/* check if current cell is obstacle cell */
	if ((obstacle & OBSTACLE_SELF) == 1){
		/* if top cell is fluid, then the no-slip boundary has to be
		 * enforced */
		if ((obstacle & OBSTACLE_TOP) == 0){
			const FLOAT dy_t = _parameters.meshsize->getDy(i,j+1);
			const FLOAT dy   = _parameters.meshsize->getDy(i,j  );
			velocity.getVector(i,j)[0] = -
				dy_t/dy*velocity.getVector(i,j+1)[0];
		}
		/* same for bottom*/
		if ((obstacle & OBSTACLE_BOTTOM) == 0){
			const FLOAT dy_b = _parameters.meshsize->getDy(i,j-1);
			const FLOAT dy   = _parameters.meshsize->getDy(i,j  );
			velocity.getVector(i,j)[0] = -
				dy_b/dy*velocity.getVector(i,j-1)[0];
		}
		/* if right cell is fluid, then the no-slip boundary has to be
		 * enforced */
		if ((obstacle & OBSTACLE_RIGHT) == 0){
			const FLOAT dx_r = _parameters.meshsize->getDx(i+1,j);
			const FLOAT dx   = _parameters.meshsize->getDx(i  ,j);
			velocity.getVector(i,j)[1] = -
				dx_r/dx*velocity.getVector(i+1,j)[1]; 
		}
		/*same for left*/
		if ((obstacle & OBSTACLE_LEFT) == 0){
			const FLOAT dx_l = _parameters.meshsize->getDx(i-1,j);
			const FLOAT dx   = _parameters.meshsize->getDx(i  ,j);
			velocity.getVector(i,j)[1] = -
				dx_l/dx*velocity.getVector(i-1,j)[1]; 
		}

		/* set normal velocity to zero if right neighbour is not
		 * obstacle*/
		if ((obstacle & OBSTACLE_RIGHT) == 0) 
			velocity.getVector(i,j)[0] = 0.0;	

		/* set normal velocity to zero if top neighbour is not
		 * obstacle*/
		if ((obstacle & OBSTACLE_TOP) == 0) 
			velocity.getVector(i,j)[1] = 0.0;	
		
	}
}

void ObstacleStencil::apply ( FlowField & flowField, int i, int j, int k){
	const int obstacle = flowField.getFlags().getValue(i, j);
	VectorField & velocity = flowField.getVelocity();    

	/* check if current cell is obstacle cell */
	if ((obstacle & OBSTACLE_SELF) == 1){
		/* if top cell is fluid: two velocities have to be set:
		 * direction 0 and 2
		 */
		if ((obstacle & OBSTACLE_TOP) == 0){
			const FLOAT dy_t = _parameters.meshsize->getDy(i,j+1,k);
			const FLOAT dy   = _parameters.meshsize->getDy(i,j  ,k);
			velocity.getVector(i,j,k)[0] = -
				dy_t/dy*velocity.getVector(i,j+1,k)[0];
			velocity.getVector(i,j,k)[2] = -
				dy_t/dy*velocity.getVector(i,j+1,k)[2];
		}
		if ((obstacle & OBSTACLE_BOTTOM) == 0){
			const FLOAT dy_b = _parameters.meshsize->getDy(i,j-1,k);
			const FLOAT dy   = _parameters.meshsize->getDy(i,j  ,k);
			velocity.getVector(i,j,k)[0] = -
				dy_b/dy*velocity.getVector(i,j-1,k)[0];
			velocity.getVector(i,j,k)[2] = -
				dy_b/dy*velocity.getVector(i,j-1,k)[2];
		}

		/* if right cell is fluid: two velocities have to be set:
		 * direction 1 and 2
		 */
		if ((obstacle & OBSTACLE_RIGHT) == 0){
			const FLOAT dx_r = _parameters.meshsize->getDx(i+1,j,k);
			const FLOAT dx   = _parameters.meshsize->getDx(i  ,j,k);
			velocity.getVector(i,j,k)[1] = -
				dx_r/dx*velocity.getVector(i+1,j,k)[1];
			velocity.getVector(i,j,k)[2] = -
				dx_r/dx*velocity.getVector(i+1,j,k)[2];
		}
		if ((obstacle & OBSTACLE_LEFT) == 0){
			const FLOAT dx_l = _parameters.meshsize->getDx(i-1,j,k);
			const FLOAT dx   = _parameters.meshsize->getDx(i  ,j,k);
			velocity.getVector(i,j,k)[1] = -
				dx_l/dx*velocity.getVector(i-1,j,k)[1];
			velocity.getVector(i,j,k)[2] = -
				dx_l/dx*velocity.getVector(i-1,j,k)[2];
		}

		/* same for fluid cell in front*/
		if ((obstacle & OBSTACLE_BACK) == 0){
			const FLOAT dz_f = _parameters.meshsize->getDx(i,j,k+1);
			const FLOAT dz   = _parameters.meshsize->getDx(i,j,k  );
			velocity.getVector(i,j,k)[1] = -
				dz_f/dz*velocity.getVector(i,j,k+1)[1];
			velocity.getVector(i,j,k)[0] = -
				dz_f/dz*velocity.getVector(i,j,k+1)[0];
		}
		if ((obstacle & OBSTACLE_FRONT) == 0){
			const FLOAT dz_b = _parameters.meshsize->getDx(i,j,k-1);
			const FLOAT dz   = _parameters.meshsize->getDx(i,j,k  );
			velocity.getVector(i,j,k)[1] = -
				dz_b/dz*velocity.getVector(i,j,k-1)[1];
			velocity.getVector(i,j,k)[0] = -
				dz_b/dz*velocity.getVector(i,j,k-1)[0];
		}

		/** Now the normal velocities need to be set to zero to ensure
		 * no flow at interfaces between solid and fluid*/
		if ((obstacle & OBSTACLE_RIGHT) == 0)
			velocity.getVector(i,j,k)[0] = 0.0;
		if ((obstacle & OBSTACLE_TOP) == 0)
			velocity.getVector(i,j,k)[1] = 0.0;
		if ((obstacle & OBSTACLE_BACK) == 0)
			velocity.getVector(i,j,k)[2] = 0.0;
	}

}
