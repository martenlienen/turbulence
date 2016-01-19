#ifndef _NSEOF_PLOTTING_LAMBDA_READER_H_
#define _NSEOF_PLOTTING_LAMBDA_READER_H_

#include <array>
#include <functional>
#include <string>

#include "Reader.h"

namespace nseof {

namespace plotting {

template <typename T, typename FF, int n>
class LambdaReader : public Reader {
 public:
  LambdaReader(FF& flowField, std::string name,
               std::function<std::array<T, n>(FF&, int, int, int)> read)
      : flowField(flowField), name(name), read(read){};

  int getDim();
  std::string getName();
  void write(const hid_t dataset, Parameters& params, nseof::hdf5::HDF5& hdf5);

 private:
  FF& flowField;
  std::string name;
  std::function<std::array<T, n>(FF&, int, int, int)> read;

  static const int LOW_OFFSET = 2;
};

template <typename T, typename FF, int n>
int LambdaReader<T, FF, n>::getDim() {
  return n;
}

template <typename T, typename FF, int n>
std::string LambdaReader<T, FF, n>::getName() {
  return this->name;
}

template <typename T, typename FF, int n>
void LambdaReader<T, FF, n>::write(const hid_t dataset, Parameters& params,
                                   nseof::hdf5::HDF5& hdf5) {
  GeometricParameters& gp = params.geometry;
  ParallelParameters& pp = params.parallel;

  int cellsX = pp.localSize[0];
  int cellsY = pp.localSize[1];
  int cellsZ = gp.dim == 3 ? pp.localSize[2] : 1;

  int ncells = cellsX * cellsY * cellsZ;

  std::vector<std::array<T, n>> buffer(ncells);

  const int low = LambdaReader::LOW_OFFSET;
  const int lowZ = gp.dim == 3 ? low : 0;

  for (int k = lowZ, a = 0; k < cellsZ + lowZ; k++) {
    for (int j = low; j < cellsY + low; j++) {
      for (int i = low; i < cellsX + low; i++, a++) {
        buffer[a] = this->read(this->flowField, i, j, k);
      }
    }
  }

  hdf5.write(dataset, buffer.data(), H5T_NATIVE_DOUBLE);
}
}
}

#endif  // _NSEOF_PLOTTING_LAMBDA_READER_H_
