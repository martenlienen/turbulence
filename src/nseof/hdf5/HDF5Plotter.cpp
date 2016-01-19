#include <array>
#include <string>
#include <sstream>
#include <vector>

#include "HDF5Plotter.h"

namespace nseof {

namespace hdf5 {

HDF5Plotter::HDF5Plotter(Parameters& params, int rank, int nranks)
    : params(params), rank(rank), hdf5(new HDF5(params.hdf5.file)) {
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

  this->xdmf->write(file.str());
}

void HDF5Plotter::plotFlowField(int timestep, const FlowField& flowField) {
  auto& readers = flowField.getReaders();

  if (this->rank == 0) {
    this->xdmf->addTimestep(timestep, readers);
  }

  std::ostringstream group;
  group << "/Data/Time-" << timestep;
  this->hdf5->createGroup(group.str());

  group << "/Rank-" << this->rank;
  this->hdf5->createGroup(group.str());

  for (auto& reader : readers) {
    std::ostringstream location;
    location << group.str() << "/" << reader->getName();

    reader->write(location.str(), this->params, *this->hdf5);
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

  this->hdf5->writeGeometry(this->rank, buffer);
}
}
}
