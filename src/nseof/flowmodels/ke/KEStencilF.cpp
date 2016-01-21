#include <algorithm>
#include <math.h>

#include "../../stencils/StencilFunctions.h"

#include "KEStencilF.h"

namespace nseof {

namespace flowmodels {

namespace ke {

KEStencilF::KEStencilF(const Parameters& parameters)
    : FieldStencil<FlowField>(parameters), cmu(parameters.kEpsilon.cmu) {
  // different wall models: chosen during initalising
  // model == 0 --> Lam and Bremhorst (Griebel)
  // model == 1 --> Chien
  // model == 2 --> Myong and Kasagi
  // model == 3 --> Fan (Paper Pennsylvania State University)
  switch (_parameters.kEpsilon.model) {
    case 1:
      calc = [this]() {
        // Lam and Bremhorst
        FLOAT Rt = wip.Rt;
        FLOAT Rd = wip.Rd;

        wop.fmu = pow(1 - exp(-0.0165 * Rd), 2) * (1 + 20.5 / max(Rt, err));
        wop.f1 = 1 + pow(0.05 / max(wop.fmu, err), 3);
        wop.f2 = 1 - exp(-Rt * Rt);
        // wop.D   = 0.0;
        // wop.E   = 0.0;
      };
      break;
    case 2:
      calc = [this]() {
        FLOAT yPlus = wip.yplus;
        FLOAT Rt = wip.Rt;
        FLOAT Re = wip.Re;
        FLOAT delta = wip.delta;
        FLOAT tke = wip.tke;
        FLOAT epsilon = wip.epsilon;

        wop.fmu = 1 - exp(-0.0115 * yPlus);
        // wop.f1  = 1.0;
        wop.f2 = 1 - (2 / 9) * exp(-Rt * Rt / 36);
        wop.D = 2 * tke / max((Re * delta * delta), err);
        wop.E =
            -2 * epsilon * exp(-0.5 * yPlus) / max((Re * delta * delta), err);
      };
      break;
    case 3:
      calc = [this]() {
        FLOAT Rt = wip.Rt;
        FLOAT yplus = wip.yplus;

        wop.fmu =
            (1 + 3.45 * pow(max(Rt, err), -0.5)) * (1 - exp(-yplus / 5.0));
        // wop.f1  = 1.0;
        wop.f2 = (1 - (2 / 9) * exp(-Rt * Rt / 36)) *
                 pow((1 - exp(-yplus / 5.0)), 2);
        // wop.D   = 0.0;
        // wop.E   = 0.0;
      };
      break;
    case 4:
      calc = [this]() {
        FLOAT Rd = wip.Rd;
        FLOAT Rt = wip.Rt;

        FLOAT fw = 1.0 - exp(-sqrt(Rd) / 2.3 +
                             (sqrt(Rd) / 2.3 - Rd / 8.89) *
                                 pow(1.0 - exp(-Rd / 20.0), 3.0));
        wop.fmu = 0.4 * fw / max(sqrt(fabs(Rt)), err) +
                  (1.0 - 0.4 * fw / max(sqrt(fabs(Rt)), err)) *
                      pow(1.0 - exp(-Rd / 42.63), 3.0);
        wop.f1 = 1.0;
        wop.f2 = (1.0 - 0.4 * exp(-Rt * Rt / 36.0) / 1.8) * fw * fw;
        wop.D = 0.0;
        wop.E = 0.0;
      };
      break;
    default:
      // do nothing
      calc = []() {};
  }

  wip.Re = parameters.flow.Re;
}

void KEStencilF::apply(FlowField& flowField, int i, int j) {
  apply(flowField, i, j, 0);
}

void KEStencilF::apply(FlowField& flowField, int i, int j, int k) {
  // turbulent kinetic energy
  wip.tke = flowField.getTke(i, j, k);
  wip.tkes = (wip.tke > 0) - (wip.tke < 0);
  wip.tke = fabs(wip.tke);

  // dissipation rate
  wip.epsilon = flowField.getEpsilon(i, j, k);
  wip.epsilons = (wip.epsilon > 0) - (wip.epsilon < 0);
  wip.epsilon = fabs(wip.epsilon);

  // closest wall distance
  wip.delta = flowField.getH(i, j, k);

  // turbulent Reynolds numbers
  wip.Rt = wip.tke * wip.tke * wip.Re / max(wip.epsilon, err);
  wip.Rd = sqrt(wip.tke) * wip.delta * wip.Re;

  // dimensionless wall distance
  FLOAT uTau = pow(cmu, 0.25) * sqrt(wip.tke);
  wip.yplus = uTau * wip.delta * wip.Re;

  // calculate factors
  calc();

  // limitation for diffusion and reaction terms
//   FLOAT gamma = flowField.getEpsilon(i,j,k)/max(1e-8, flowField.getTke(i,j,k));
//   wop.f2 *= (wop.f2*gamma> 0.0) ? 1.0 : 0.0;
//   wop.f3  = (       gamma> 0.0) ? 1.0 : 0.0; 

  // write factors
  flowField.getFmu(i, j, k) = wop.fmu;
  flowField.getF1(i, j, k) = wop.f1;
  flowField.getF2(i, j, k) = wop.f2;
  flowField.getF3(i, j, k) = wop.f3;
  flowField.getD(i, j, k) = wop.D;
  flowField.getE(i, j, k) = wop.E;
}
}
}
}
