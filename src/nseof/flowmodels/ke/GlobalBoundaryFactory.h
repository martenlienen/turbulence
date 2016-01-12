#ifndef _NSEOF_FLOWMODELS_KE_GLOBAL_BOUNDARY_FACTORY_H_
#define _NSEOF_FLOWMODELS_KE_GLOBAL_BOUNDARY_FACTORY_H_

#include <string>

#include "../../GlobalBoundaryFactory.h"
#include "../../Parameters.h"
#include "../../Iterators.h"

#include "../../stencils/MovingWallStencils.h"
#include "../../stencils/PeriodicBoundaryStencils.h"
#include "../../stencils/NeumannBoundaryStencils.h"
#include "../../stencils/BFInputStencils.h"

#include "FlowField.h"

namespace nseof {

namespace flowmodels {

namespace ke {

class GlobalBoundaryFactory : public nseof::GlobalBoundaryFactory {
 private:
  BoundaryStencil<FlowField>* _keStencils[6];

 public:
  GlobalBoundaryFactory(Parameters& parameters);

  ~GlobalBoundaryFactory();

  GlobalBoundaryIterator<FlowField> getGlobalBoundaryIterator(
      FlowField& flowField);
};
}
}
}

#endif
