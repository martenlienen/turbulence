#ifndef _NSEOF_FLOWMODELS_ALGEBRAIC_FLOWFIELD_H_
#define _NSEOF_FLOWMODELS_ALGEBRAIC_FLOWFIELD_H_

#include "../turbulent/FlowField.h"

namespace nseof {

namespace flowmodels {

namespace algebraic {

/**
 * Flow field for an algebraic turbulence model
 */
class FlowField : public nseof::flowmodels::turbulent::FlowField {
 public:
  FlowField(const Parameters& parameters)
      : nseof::flowmodels::turbulent::FlowField(parameters) {}
};
}
}
}

#endif
