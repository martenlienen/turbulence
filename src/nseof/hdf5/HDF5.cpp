// Ensure that mpi.h is included before everything else (otherwise MAC-Cluster
// complains)
#include "HDF5.h"

#include <sstream>

namespace nseof {

namespace hdf5 {

HDF5::HDF5(std::string path) { this->open(path); }
HDF5::~HDF5() { this->close(); }

void HDF5::write(const hid_t dataset, const void* buffer, hid_t type) {
  hid_t plist = H5Pcreate(H5P_DATASET_XFER);
  H5Pset_dxpl_mpio(plist, H5FD_MPIO_INDEPENDENT);

  H5Dwrite(dataset, type, H5S_ALL, H5S_ALL, plist, buffer);

  H5Pclose(plist);
}

void HDF5::createGroup(std::string location) {
  hid_t group = H5Gcreate2(this->file, location.c_str(), H5P_DEFAULT,
                           H5P_DEFAULT, H5P_DEFAULT);
  H5Gclose(group);
}

const hid_t& HDF5::getFile() { return this->file; }
const hid_t& HDF5::getGeometry() { return this->geometry; }

void HDF5::open(std::string path) {
  hid_t plist = H5Pcreate(H5P_FILE_ACCESS);
  H5Pset_fapl_mpio(plist, MPI_COMM_WORLD, MPI_INFO_NULL);

  this->file = H5Fcreate(path.c_str(), H5F_ACC_TRUNC, H5P_DEFAULT, plist);
  H5Pclose(plist);

  this->geometry = H5Gcreate2(this->file, "Geometries", H5P_DEFAULT,
                              H5P_DEFAULT, H5P_DEFAULT);
  this->data =
      H5Gcreate2(this->file, "Data", H5P_DEFAULT, H5P_DEFAULT, H5P_DEFAULT);
}

void HDF5::close() {
  H5Gclose(this->data);
  H5Gclose(this->geometry);
  H5Fclose(this->file);
}
}
}
