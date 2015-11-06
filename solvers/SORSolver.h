#ifndef _SOLVER_H_
#define _SOLVER_H_

// Yhis moduke solves the diffusion equation for the pressure. Right now it's very simple.

#include "../FlowField.h"
#include "../Definitions.h"
#include "../Parameters.h"
#include "../LinearSolver.h"

class SORSolver : public LinearSolver {
    public:
        SORSolver(FlowField & flowField, const Parameters & parameters);
        void solve ();
};

#endif
