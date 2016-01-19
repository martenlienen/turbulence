#ifndef _NSEOF_HDF5_HDF5_H_
#define _NSEOF_HDF5_HDF5_H_

#include <array>
#include <string>
#include <vector>

#include <hdf5.h>

namespace nseof {

namespace hdf5 {

class HDF5 {
 public:
  HDF5(std::string path);
  ~HDF5();

  void writeGeometry(int rank, const std::vector<std::array<float, 3>>& points);

  void write(std::string location, hsize_t n, hsize_t dim, const void* buffer,
             hid_t type);

  void createGroup(std::string location);

 private:
  void open(std::string path);
  void close();

  hid_t file;
  hid_t geometry;
  hid_t data;
};
}
}

#endif  // _NSEOF_HDF5_HDF5_H_
