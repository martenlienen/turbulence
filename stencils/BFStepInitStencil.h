#ifndef _BF_STEP_INIT_STENCIL_H_
#define _BF_STEP_INIT_STENCIL_H_

#include "../Stencil.h"
#include "../FlowField.h"
#include "../Definitions.h"

/** initialized the backward facing step scenario, i.e. sets the flag field.
 *
 */
class BFStepInitStencil : public FieldStencil<FlowField> {
    public:

        BFStepInitStencil (const Parameters & parameters);

        void apply(FlowField & flowField, int i, int j);
        void apply(FlowField & flowField, int i, int j, int k);

    private:

        const FLOAT xLimit;   //! size of step in x-direction
        const FLOAT yLimit;   //! Same as for x
};

#endif
