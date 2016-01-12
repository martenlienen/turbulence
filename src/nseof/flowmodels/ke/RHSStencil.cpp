#include <functional>

#include "../../stencils/StencilFunctions.h"

#include "FlowField.h"
#include "RHSStencil.h"

namespace nseof {

namespace flowmodels {

namespace ke {

RHSStencil::RHSStencil(const Parameters& parameters)
    : FieldStencil<FlowField>(parameters) {}

void RHSStencil::apply(FlowField& flowField, int i, int j) {
  //  loadLocalNu2D(_parameters, flowField, _localNu, i, j);
  loadLocalVelocity2D(flowField, _localVelocity, i, j);
  loadLocalMeshsize2D(_parameters, _localMeshsize, i, j);

  // load TKE
  loadLocal2D([&flowField](FLOAT* local, int ii, int jj) mutable {
    *(local + 0) = flowField.getTke(ii, jj);
  }, _localTKE, i, j);

  // load epsilon
  loadLocal2D([&flowField](FLOAT* local, int ii, int jj) mutable {
    *(local + 0) = flowField.getEpsilon(ii, jj);
  }, _localEpsilon, i, j);

  // load Fmu*Nut
  loadLocal2D([&flowField](FLOAT* local, int ii, int jj) mutable {
    *(local + 0) = flowField.getFmu(ii, jj) * flowField.getNu(ii, jj);
    //    *(local+0) = 0;
  }, _localFmuNut, i, j);

  const FLOAT sijsij = computeSijSij2D(_localVelocity, _localMeshsize);
  const FLOAT f1 = flowField.getF1(i, j);
  const FLOAT f2 = flowField.getF2(i, j);
  const FLOAT f3 = flowField.getF3(i, j);
  const FLOAT D = flowField.getD(i, j);
  const FLOAT E = flowField.getE(i, j);
  const FLOAT nut = flowField.getNu(i, j);

  flowField.getsijsij(i, j) = sijsij;

  // load nu + nut/sigmaK
  loadLocal2D([&flowField, this](FLOAT* local, int ii, int jj) mutable {
    *(local + 0) = this->_parameters.flow.visc +
                   flowField.getNu(ii, jj) / this->_parameters.kEpsilon.sigmaK;
  }, _localNu, i, j);

  flowField.getRHSTke(i, j) =
      computeRHStke(_parameters, _localVelocity, _localMeshsize, _localNu,
                    _localTKE, _localEpsilon, nut, sijsij, f3, D);

  // load nu + nut/sigmaE
  loadLocal2D([&flowField, this](FLOAT* local, int ii, int jj) mutable {
    *(local + 0) = this->_parameters.flow.visc +
                   flowField.getNu(ii, jj) / this->_parameters.kEpsilon.sigmaE;
  }, _localNu, i, j);

  flowField.getRHSEpsilon(i, j) =
      computeRHSepsilon(_parameters, _localVelocity, _localMeshsize, _localNu,
                        _localTKE, _localEpsilon, nut, sijsij, f1, f2, E);
}

void RHSStencil::apply(FlowField& flowField, int i, int j, int k) {
  // TODO: 3D
}
}
}
}
