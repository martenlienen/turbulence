#ifndef _INITTAYLORGREENFLOWFIELDSTENCIL_H_
#define _INITTAYLORGREENFLOWFIELDSTENCIL_H_

#include "../Definitions.h"
#include "../Parameters.h"
#include "../Stencil.h"
#include "../FlowField.h"
#include <cmath>

/** 
 *  Stencil for initializing the taylor-green vortex flow in a periodic domain.
 *  @author Philipp Neumann
 */
class InitTaylorGreenFlowFieldStencil : public FieldStencil<FlowField> {

    private:
      const FLOAT _2pi;
      const FLOAT * const  _domainSize;

      FLOAT * initializeDomainSize(const Parameters &parameters) const {
        FLOAT *domainSize = new FLOAT[3];
        if (domainSize == NULL){handleError(1,"InitTaylorGreenFlowFieldStencil: domainSize==NULL");}
        domainSize[0] = parameters.geometry.lengthX;
        domainSize[1] = parameters.geometry.lengthY;
        domainSize[2] = parameters.geometry.lengthZ;
        return domainSize;
      }

      /** from the local grid coordinates i,j,k, computes the global coordinates of the current cell and initializes
       *  the velocity field correspondingly
       */
      void computeGlobalCoordinates(FLOAT *coords,int i, int j,int k=0) const {
        coords[0] = _parameters.meshsize->getPosX(i,j,k) + 0.5*_parameters.meshsize->getDx(i,j,k);
        coords[1] = _parameters.meshsize->getPosY(i,j,k) + 0.5*_parameters.meshsize->getDy(i,j,k);
        if (_parameters.geometry.dim == 3){
          coords[2] = _parameters.meshsize->getPosZ(i,j,k) + 0.5*_parameters.meshsize->getDz(i,j,k);
        }
        // Output for debugging
        // std::cout << "RankCorner " << _parameters.parallel.rank << ": " << _parameters.parallel.firstCorner[0] << "," << _parameters.parallel.firstCorner[1] << "," << _parameters.parallel.firstCorner[2] << std::endl;
        // std::cout << "Rank " << _parameters.parallel.rank << ": " << coords[0] << "," << coords[1] << "," << coords[2] << std::endl;
      }

    public:

        /** Constructor
         *
         */
        InitTaylorGreenFlowFieldStencil ( const Parameters & parameters ): FieldStencil<FlowField>(parameters),
          _2pi(2.0*3.141592653589793238), _domainSize(initializeDomainSize(parameters))
        {}

        virtual ~InitTaylorGreenFlowFieldStencil(){
          if (_domainSize != NULL){delete [] _domainSize;}
        }

        /** 2D operation for one position
         *
         * @param flowField State of the flow field
         * @param i Position in the x direction
         * @param j Position in the y direction
         */
        void apply ( FlowField & flowField, int i, int j ){
          FLOAT coords[3]={0.0,0.0,0.0};
          FLOAT * const velocity = flowField.getVelocity().getVector(i,j);
          computeGlobalCoordinates(coords,i,j);
          // initialize velocities
          velocity[0] = sin(_2pi*(coords[0]+0.5*_parameters.meshsize->getDx(i,j))/_domainSize[0])*
                        sin(_2pi*coords[1]/_domainSize[1]);
          velocity[1] = cos(_2pi*coords[0]/_domainSize[0])*
                        cos(_2pi*(coords[1]+0.5*_parameters.meshsize->getDy(i,j))/_domainSize[1]);
        }

        /** 3D operation for one position
         *
         * @param flowField State of the flow field
         * @param i Position in the x direction
         * @param j Position in the y direction
         * @param k Position in the z direction
         */
        void apply ( FlowField & flowField, int i, int j, int k ){
          FLOAT coords[3]={0.0,0.0,0.0};
          FLOAT * const velocity = flowField.getVelocity().getVector(i,j,k);
          computeGlobalCoordinates(coords,i,j,k);
          // initialize velocities
          velocity[0] = cos(_2pi*(coords[0]+0.5*_parameters.meshsize->getDx(i,j,k))/_domainSize[0])*
                        sin(_2pi*coords[1]/_domainSize[1])*
                        sin(_2pi*coords[2]/_domainSize[2]);
          velocity[1] = sin(_2pi*coords[0]/_domainSize[0])*
                        cos(_2pi*(coords[1]+0.5*_parameters.meshsize->getDy(i,j,k))/_domainSize[1])*
                        sin(_2pi*coords[2]/_domainSize[2]);
          velocity[2] = sin(_2pi*coords[0]/_domainSize[0])*
                        sin(_2pi*coords[1]/_domainSize[1])*
                        cos(_2pi*(coords[2]+0.5*_parameters.meshsize->getDz(i,j,k))/_domainSize[2]);
        }
};

#endif
