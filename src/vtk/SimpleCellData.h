#ifndef _VTK_SIMPLE_CELL_DATA_H
#define _VTK_SIMPLE_CELL_DATA_H

#include <vector>

#include "CellData.h"

namespace vtk {

template <typename T>
class SimpleCellData : public CellData {
 public:
  SimpleCellData(std::string dataName, std::vector<T> data);
  std::string str();

 private:
  std::vector<T> data;
};
}

// Implementation

#include <ios>
#include <sstream>

namespace vtk {

template <typename T>
SimpleCellData<T>::SimpleCellData(std::string dataName, std::vector<T> data)
    : CellData(dataName), data(data) {}

template <>
std::string SimpleCellData<double>::str() {
  std::ostringstream stream;
  stream << std::fixed;

  stream << "SCALARS " << this->getDataName() << " float 1\n";
  stream << "LOOKUP_TABLE default\n";

  for (auto& scalar : this->data) {
    stream << scalar << "\n";
  }

  return stream.str();
}

template <>
std::string SimpleCellData<std::vector<double>>::str() {
  std::ostringstream stream;
  stream << std::fixed;

  stream << "VECTORS " << this->getDataName() << " float\n";

  for (auto& vector : this->data) {
    bool first = true;

    for (auto& scalar : vector) {
      if (first) {
        first = false;
      } else {
        stream << " ";
      }

      stream << scalar;
    }

    stream << "\n";
  }

  return stream.str();
}
}

#endif  // _VTK_SIMPLE_CELL_DATA_H
