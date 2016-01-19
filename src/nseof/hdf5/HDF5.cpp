#include <sstream>

#include "HDF5.h"

namespace nseof {

namespace hdf5 {

HDF5::HDF5(std::string path) { this->open(path); }
HDF5::~HDF5() { this->close(); }

void HDF5::writeGeometry(int rank,
                         const std::vector<std::array<float, 3>>& points) {
  const hsize_t dimensions[2] = {points.size(), 3};
  hid_t dataspace = H5Screate_simple(2, dimensions, NULL);

  std::ostringstream name;
  name << "Rank-" << rank;
  hid_t dataset =
      H5Dcreate(this->geometry, name.str().c_str(), H5T_NATIVE_FLOAT, dataspace,
                H5P_DEFAULT, H5P_DEFAULT, H5P_DEFAULT);

  H5Dwrite(dataset, H5T_NATIVE_FLOAT, H5S_ALL, H5S_ALL, H5P_DEFAULT,
           points.data());

  H5Dclose(dataset);
  H5Sclose(dataspace);
}

void HDF5::write(std::string location, hsize_t n, hsize_t dim,
                 const void* buffer, hid_t type) {
  const hsize_t dimensions[2] = {n, dim};
  hid_t dataspace = H5Screate_simple(2, dimensions, NULL);

  hid_t dataset = H5Dcreate(this->data, location.c_str(), type, dataspace,
                            H5P_DEFAULT, H5P_DEFAULT, H5P_DEFAULT);

  H5Dwrite(dataset, type, H5S_ALL, H5S_ALL, H5P_DEFAULT, buffer);

  H5Dclose(dataset);
  H5Sclose(dataspace);
}

void HDF5::createGroup(std::string location) {
  hid_t group = H5Gcreate2(this->file, location.c_str(), H5P_DEFAULT,
                           H5P_DEFAULT, H5P_DEFAULT);
  H5Gclose(group);
}

void HDF5::open(std::string path) {
  this->file = H5Fcreate(path.c_str(), H5F_ACC_EXCL, H5P_DEFAULT, H5P_DEFAULT);
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
