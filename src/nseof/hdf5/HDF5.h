#ifndef _NSEOF_HDF5_HDF5_H_
#define _NSEOF_HDF5_HDF5_H_

// Ensure that mpi.h is included before everything else (otherwise MAC-Cluster
// complains)
#include <hdf5.h>

#include <array>
#include <string>
#include <vector>

namespace nseof {

namespace hdf5 {

class HDF5 {
 public:
  HDF5(std::string path);
  ~HDF5();

  void write(const hid_t dataset, const void* buffer, hid_t type);

  void createGroup(std::string location);

  const hid_t& getFile();
  const hid_t& getGeometry();

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
