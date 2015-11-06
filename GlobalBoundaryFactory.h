#ifndef _FACTORIES_H_
#define _FACTORIES_H_

#include <string>
#include "Parameters.h"
#include "Iterators.h"
#include "stencils/MovingWallStencils.h"
#include "stencils/PeriodicBoundaryStencils.h"
#include "stencils/NeumannBoundaryStencils.h"
#include "stencils/BFInputStencils.h"
#include "FlowField.h"

/** Class that returns instances of the global boundary iterator. It also contains the stencils.
 * Right now, it works only with Dirichlet and periodic boundary conditions
 */
class GlobalBoundaryFactory{

    private:
        // List of all the stencils
        BoundaryStencil<FlowField> * _velocityStencils[6];    //! A stencil for each face
        BoundaryStencil<FlowField> * _FGHStencils[6];    //! A stencil for each face
        BoundaryStencil<FlowField> * _moving[2];        //! Pointers to the moving wall stencils, if any
        BoundaryStencil<FlowField> * _periodic[2];      //! Pointers to the periodic stencils, if any
        BoundaryStencil<FlowField> * _outflow[2];     //! Pointers for the outflow conditions
        BoundaryStencil<FlowField> * _channelInput[2];    //! For the velocity input
        const Parameters & _parameters;         //! Reference to the parameters

    public:

        /** Constructor. Will initialize all the references to the boundaries and create the
         * stencils.
         * @param parameters Parameters of the problem
         */
        GlobalBoundaryFactory(Parameters & parameters);

        /** Destructor. Necessary to remove the stencils at the end
         */
        ~GlobalBoundaryFactory();

        /** Returns an instance of the global boundary iterator for velocities.
         * @param flowField Flow field information
         */
        GlobalBoundaryIterator<FlowField> getGlobalBoundaryVelocityIterator(FlowField & flowField);

        /** Returns an instance of the global boundary iterator for FGH.
         * @param flowField Flow field information
         */
        GlobalBoundaryIterator<FlowField> getGlobalBoundaryFGHIterator(FlowField & flowField);
};

#endif
