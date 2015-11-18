#include "File.h"

#include <iostream>
#include <fstream>
#include <sstream>

namespace vtk {

File::File(Dataset dataset, std::vector<std::unique_ptr<CellData>> cellData)
    : dataset(dataset), cellData(std::move(cellData)) {}

bool File::write(std::string prefix, int timestep) {
  std::ostringstream path;

  path << prefix << "." << timestep << ".vtk";

  return this->write(path.str());
}

bool File::write(std::string path) {
  std::fstream file;
  file.open(path.c_str(), std::ios::out);

  if (file.fail()) {
    return false;
  }

  file << this->header();
  file << this->dataset.str();

  if (!this->cellData.empty()) {
    file << this->cellDataHeader();

    for (auto& data : this->cellData) {
      file << data->str();
    }
  }

  return true;
}

std::string File::header() {
  std::ostringstream header;

  header << "# vtk DataFile Version 2.0\n";
  header << "NS-EOF\n";
  header << "ASCII\n";

  return header.str();
}

std::string File::cellDataHeader() {
  std::ostringstream header;

  header << "CELL_DATA " << this->dataset.numCells() << "\n";

  return header.str();
}
}
