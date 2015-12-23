#ifndef _VTK_FILE_H
#define _VTK_FILE_H

#include <memory>
#include <string>
#include <vector>

#include "Dataset.h"
#include "CellData.h"

namespace vtk {

class File {
 public:
  File(Dataset dataset, std::vector<std::unique_ptr<CellData>> cellData);

  bool write(std::string prefix, int rank, int timestep);
  bool write(std::string prefix, int timestep);
  bool write(std::string path);

 private:
  Dataset dataset;

  std::vector<std::unique_ptr<CellData>> cellData;

  std::string header();
  std::string cellDataHeader();
};
}

#endif  // _VTK_FILE_H
