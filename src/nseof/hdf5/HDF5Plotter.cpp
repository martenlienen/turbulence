#include <array>
#include <string>
#include <sstream>
#include <vector>

#include "HDF5Plotter.h"

namespace nseof {

namespace hdf5 {

HDF5Plotter::HDF5Plotter(Parameters& params, int rank, int nranks)
    : params(params),
      rank(rank),
      nranks(nranks),
      hdf5(new HDF5(params.hdf5.file)) {
  if (rank == 0) {
    this->xdmf = XDMF::fromParameters(params, params.hdf5.file, nranks);
  }

  this->plotGeometry(params);
}

HDF5Plotter::~HDF5Plotter() {
  std::ostringstream file;
  std::string h5file = this->params.hdf5.file;

  // If h5file ends with .h5
  size_t pos = h5file.rfind(".h5");
  if (pos == h5file.size() - 3) {
    file << h5file.substr(0, pos) << ".xmf";
  } else {
    file << h5file << ".xmf";
  }

  if (this->rank == 0) {
    this->xdmf->write(file.str());
  }
}

void HDF5Plotter::plotFlowField(int timestep, const FlowField& flowField) {
  auto& readers = flowField.getReaders();

  if (this->rank == 0) {
    this->xdmf->addTimestep(timestep, readers);
  }

  GeometricParameters& gp = this->params.geometry;
  ParallelParameters& pp = this->params.parallel;

  int cellsX = pp.localSize[0];
  int cellsY = pp.localSize[1];
  int cellsZ = gp.dim == 3 ? pp.localSize[2] : 1;
  int ncells = cellsX * cellsY * cellsZ;

  std::ostringstream group;
  group << "/Data/Time-" << timestep;
  this->hdf5->createGroup(group.str());

  std::vector<std::vector<hid_t>> datasets(this->nranks,
                                           std::vector<hid_t>(readers.size()));

  for (int i = 0; i < this->nranks; i++) {
    std::ostringstream rankGroup;
    rankGroup << group.str() << "/Rank-" << i;
    this->hdf5->createGroup(rankGroup.str());

    for (size_t j = 0; j < readers.size(); j++) {
      std::ostringstream location;
      location << rankGroup.str() << "/" << readers[j]->getName();

      const hsize_t dimensions[2] = {(hsize_t)ncells,
                                     (hsize_t)readers[j]->getDim()};
      hid_t dataspace = H5Screate_simple(2, dimensions, NULL);

      hid_t dataset = H5Dcreate(this->hdf5->getFile(), location.str().c_str(),
                                readers[j]->getHDF5NativeType(), dataspace,
                                H5P_DEFAULT, H5P_DEFAULT, H5P_DEFAULT);

      H5Sclose(dataspace);

      datasets[i][j] = dataset;
    }
  }

  for (size_t i = 0; i < readers.size(); i++) {
    readers[i]->write(datasets[this->rank][i], this->params, *this->hdf5);
  }

  for (int i = 0; i < this->nranks; i++) {
    for (size_t j = 0; j < readers.size(); j++) {
      H5Dclose(datasets[i][j]);
    }
  }
}

void HDF5Plotter::plotGeometry(Parameters& params) {
  GeometricParameters& gp = params.geometry;
  ParallelParameters& pp = params.parallel;
  Meshsize* mesh = params.meshsize;

  int cellsX = pp.localSize[0];
  int cellsY = pp.localSize[1];
  int cellsZ = gp.dim == 3 ? pp.localSize[2] : 1;

  int pointsX = cellsX + 1;
  int pointsY = cellsY + 1;
  int pointsZ = gp.dim == 3 ? cellsZ + 1 : 1;

  int npoints = pointsX * pointsY * pointsZ;

  std::vector<std::array<float, 3>> buffer(npoints);

  const int low = HDF5Plotter::LOW_OFFSET;
  const int lowZ = gp.dim == 3 ? low : 0;

  for (int z = lowZ, i = 0; z < pointsZ + lowZ; z++) {
    for (int y = low; y < pointsY + low; y++) {
      for (int x = low; x < pointsX + low; x++, i++) {
        buffer[i][0] = mesh->getPosX(x, y, z);
        buffer[i][1] = mesh->getPosY(x, y, z);
        buffer[i][2] = mesh->getPosZ(x, y, z);
      }
    }
  }

  std::vector<hid_t> datasets(this->nranks);

  for (int i = 0; i < this->nranks; i++) {
    const hsize_t dimensions[2] = {(hsize_t)npoints, 3};
    hid_t dataspace = H5Screate_simple(2, dimensions, NULL);

    std::ostringstream name;
    name << "Rank-" << i;
    hid_t dataset = H5Dcreate(this->hdf5->getGeometry(), name.str().c_str(),
                              H5T_NATIVE_FLOAT, dataspace, H5P_DEFAULT,
                              H5P_DEFAULT, H5P_DEFAULT);

    H5Sclose(dataspace);

    datasets[i] = dataset;
  }

  this->hdf5->write(datasets[this->rank], buffer.data(), H5T_NATIVE_FLOAT);

  for (auto& dataset : datasets) {
    H5Dclose(dataset);
  }
}
}
}
