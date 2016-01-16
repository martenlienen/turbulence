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
  //  loadLocalNu3D(_parameters, flowField, _localNu, i, j, k);
  loadLocalVelocity3D(flowField, _localVelocity, i, j, k);
  loadLocalMeshsize3D(_parameters, _localMeshsize, i, j, k);

  // load TKE
  loadLocal3D([&flowField](FLOAT* local, int ii, int jj, int kk) mutable {
    *(local + 0) = flowField.getTke(ii, jj, kk);
  }, _localTKE, i, j, k);

  // load epsilon
  loadLocal3D([&flowField](FLOAT* local, int ii, int jj, int kk) mutable {
    *(local + 0) = flowField.getEpsilon(ii, jj, kk);
  }, _localEpsilon, i, j, k);

  // load Fmu*Nut
  loadLocal3D([&flowField](FLOAT* local, int ii, int jj, int kk) mutable {
    *(local + 0) = flowField.getFmu(ii, jj, kk) * flowField.getNu(ii, jj, kk);
    //    *(local+0) = 0;
  }, _localFmuNut, i, j, k);

  const FLOAT sijsij = computeSijSij3D(_localVelocity, _localMeshsize);
  const FLOAT f1 = flowField.getF1(i, j, k);
  const FLOAT f2 = flowField.getF2(i, j, k);
  const FLOAT f3 = flowField.getF3(i, j, k);
  const FLOAT D = flowField.getD(i, j, k);
  const FLOAT E = flowField.getE(i, j, k);
  const FLOAT nut = flowField.getNu(i, j, k);

  flowField.getsijsij(i, j, k) = sijsij;

  // load nu + nut/sigmaK
  loadLocal3D([&flowField, this](FLOAT* local, int ii, int jj, int kk) mutable {
    *(local + 0) = this->_parameters.flow.visc +
                   flowField.getNu(ii, jj, kk) / this->_parameters.kEpsilon.sigmaK;
  }, _localNu, i, j, k);

  flowField.getRHSTke(i, j, k) =
      computeRHStke(_parameters, _localVelocity, _localMeshsize, _localNu,
                    _localTKE, _localEpsilon, nut, sijsij, f3, D);

  // load nu + nut/sigmaE
  loadLocal3D([&flowField, this](FLOAT* local, int ii, int jj, int kk) mutable {
    *(local + 0) = this->_parameters.flow.visc +
                   flowField.getNu(ii, jj, kk) / this->_parameters.kEpsilon.sigmaE;
  }, _localNu, i, j, k);

  flowField.getRHSEpsilon(i, j, k) =
      computeRHSepsilon(_parameters, _localVelocity, _localMeshsize, _localNu,
                        _localTKE, _localEpsilon, nut, sijsij, f1, f2, E);
}
}
}
}
