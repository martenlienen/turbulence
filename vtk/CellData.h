#ifndef _VTK_CELL_DATA_H
#define _VTK_CELL_DATA_H

#include <string>

namespace vtk {

class CellData {
 public:
  CellData(std::string dataName);
  virtual ~CellData(){};

  virtual std::string str() = 0;

  std::string getDataName();

 private:
  std::string dataName;
};
}

#endif  // _VTK_CELL_DATA_H
