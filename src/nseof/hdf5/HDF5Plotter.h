#ifndef _NSEOF_HDF5_HDF5PLOTTER_H_
#define _NSEOF_HDF5_HDF5PLOTTER_H_

#include <memory>
#include <string>

#include "../Parameters.h"
#include "../FlowField.h"

#include "HDF5.h"
#include "XDMF.h"

namespace nseof {

namespace hdf5 {

class HDF5Plotter {
 public:
  static const int LOW_OFFSET = 2;
  static const int HIGH_OFFSET = 1;

  HDF5Plotter(Parameters& params, int rank, int nranks);
  ~HDF5Plotter();

  void plotFlowField(int timestep, const FlowField& flowField);

 private:
  Parameters& params;
  int rank;
  int nranks;
  std::unique_ptr<HDF5> hdf5;
  std::unique_ptr<XDMF> xdmf;

  void plotGeometry(Parameters& params);
};
}
}

#endif  // _NSEOF_HDF5_HDF5PLOTTER_H_
