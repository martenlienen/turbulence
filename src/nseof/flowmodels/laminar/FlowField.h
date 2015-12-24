#ifndef _NSEOF_FLOWMODELS_LAMINAR_FLOW_FIELD_H_
#define _NSEOF_FLOWMODELS_LAMINAR_FLOW_FIELD_H_

#include "../../DataStructures.h"
#include "../../Parameters.h"
#include "../../FlowField.h"

namespace nseof {

namespace flowmodels {

namespace laminar {

/**
 * A laminar flow field
 */
class FlowField : public nseof::FlowField {
 public:
  FlowField(const Parameters& params) : nseof::FlowField(params) {}
};
}
}
}

#endif
