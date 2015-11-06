#include "LinearSolver.h"

LinearSolver::LinearSolver(FlowField & flowField, const Parameters & parameters):
    _flowField(flowField), _parameters(parameters){}
