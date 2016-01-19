#ifndef _NSEOF_HDF5_XDMF_H_
#define _NSEOF_HDF5_XDMF_H_

#include <memory>
#include <string>
#include <memory>
#include <vector>

#include <tinyxml2.h>

#include "../Parameters.h"
#include "../plotting/Reader.h"

namespace nseof {

namespace hdf5 {

/**
 * A description of the HDF5 output
 *
 * This assumes that all subdomains are of the same size.
 */
class XDMF {
 public:
  static std::unique_ptr<XDMF> fromParameters(Parameters& params,
                                              std::string hdf, int nranks);

  /**
   * @param hdf Path to HDF file relative to XDMF file
   * @param cellsX Number of cells in subdomains in X direction
   * @param cellsY Number of cells in subdomains in Y direction
   * @param cellsZ Number of cells in subdomains in Z direction
   * @param pointsX Number of points in subdomains in X direction
   * @param pointsY Number of points in subdomains in Y direction
   * @param pointsZ Number of points in subdomains in Z direction
   * @param nranks Total number of MPI ranks
   */
  XDMF(std::string hdf, int cellsX, int cellsY, int cellsZ, int pointsX,
       int pointsY, int pointsZ, int nranks);

  void addTimestep(
      int timestep,
      const std::vector<std::unique_ptr<nseof::plotting::Reader>>& readers);

  void write(std::string filename);

 private:
  static const int GHOST_LAYERS = 3;

  std::string hdf;
  int nranks;
  int cellsX;
  int cellsY;
  int cellsZ;
  tinyxml2::XMLDocument doc;
  tinyxml2::XMLElement* timesteps;
};
}
}

#endif  // _NSEOF_HDF5_XDMF_H_
