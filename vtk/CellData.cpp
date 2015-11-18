#include "CellData.h"

namespace vtk {

CellData::CellData(std::string dataName) : dataName(dataName) {}

std::string CellData::getDataName() { return this->dataName; }
}
