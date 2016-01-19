#ifndef _NSEOF_PLOTTING_READER_H_
#define _NSEOF_PLOTTING_READER_H_

#include <string>

#include "../hdf5/HDF5.h"

namespace nseof {

namespace plotting {

class Reader {
 public:
  virtual int getDim() = 0;
  virtual std::string getName() = 0;

  /**
   * Read data from the flow field and write it into an hdf5 file
   */
  virtual void write(std::string location, Parameters& p,
                     nseof::hdf5::HDF5& hdf5) = 0;
};
}
}

#endif  // _NSEOF_PLOTTING_READER_H_
