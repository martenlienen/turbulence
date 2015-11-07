#include "VTKStencil.h"

#define GHOST_OFFSET 2

#include <sstream>
#include <fstream>
#include <algorithm>
#include <vector>

VTKStencil::VTKStencil(const Parameters& parameters)
  : FieldStencil(parameters) {
}

void VTKStencil::apply(FlowField& flowField, int i, int j) {

}

void VTKStencil::apply(FlowField& flowField, int i, int j, int k) {

}

void VTKStencil::write(FlowField& flowField, int timeStep) {
  GeometricParameters geom = this->_parameters.geometry;
  Meshsize* mesh = this->_parameters.meshsize;
  IntScalarField& flags = flowField.getFlags();

  // +3 is for the ghost layer
  int cellsX = geom.sizeX;
  int cellsY = geom.sizeY;
  int cellsZ = geom.sizeZ;
  int cells = cellsX * cellsY * cellsZ;
  int pointsX = cellsX + 1;
  int pointsY = cellsY + 1;
  int pointsZ = geom.dim == 3 ? cellsZ + 1 : 1;
  int points = pointsX * pointsY * pointsZ;
  std::vector<float> pressures(cells);
  std::vector<float> velocitiesX(cells);
  std::vector<float> velocitiesY(cells);
  std::vector<float> velocitiesZ(cells);

  for (int k = 0, cell = 0; k < cellsZ; k++) {
    for (int j = 0; j < cellsY; j++) {
      for (int i = 0; i < cellsX; i++, cell++) {
        FLOAT p;
        FLOAT v[3];
        int x = i + GHOST_OFFSET;
        int y = j + GHOST_OFFSET;
        int z = k + GHOST_OFFSET;

        if ((flags.getValue(x, y, k) & OBSTACLE_SELF) == 0) {
          if (geom.dim == 2) {
            flowField.getPressureAndVelocity(p, v, x, y);
          } else {
            flowField.getPressureAndVelocity(p, v, x, y, z);
          }
        } else {
          p = 0.0;
        }

        pressures[cell] = p;
      }
    }
  }

  for (int k = 0, cell = 0; k < cellsZ; k++) {
    for (int j = 0; j < cellsY; j++) {
      for (int i = 0; i < cellsX; i++, cell++) {
        FLOAT p;
        FLOAT v[3] = {0, 0, 0};
        int x = i + GHOST_OFFSET;
        int y = j + GHOST_OFFSET;
        int z = k + GHOST_OFFSET;

        if ((flags.getValue(x, y, k) & OBSTACLE_SELF) == 0) {
          if (geom.dim == 2) {
            flowField.getPressureAndVelocity(p, v, x, y);
          } else {
            flowField.getPressureAndVelocity(p, v, x, y, z);
          }
        }

        velocitiesX[cell] = v[0];
        velocitiesY[cell] = v[1];
        velocitiesZ[cell] = v[2];
      }
    }
  }

  std::stringstream fstream;
  fstream << _parameters.vtk.prefix << "." << timeStep << ".vtk";
  std::string filename = fstream.str();
  std::ofstream file;
  file.open(filename.c_str(), std::ios::out);

  // Print floats with fixed precision
  file << std::fixed;

  file << "# vtk DataFile Version 2.0\n";
  file << "NS-EOF\n";
  file << "ASCII\n";
  file << "DATASET STRUCTURED_GRID\n";
  file << "DIMENSIONS " << pointsX << " " << pointsY << " " << pointsZ << "\n";
  file << "POINTS " << points << " float\n";

  for (int k = GHOST_OFFSET; k < pointsZ + GHOST_OFFSET; k++) {
    for (int j = GHOST_OFFSET; j < pointsY + GHOST_OFFSET; j++) {
      for (int i = GHOST_OFFSET; i < pointsX + GHOST_OFFSET; i++) {
        FLOAT posX = mesh->getPosX(i, j, k);
        FLOAT posY = mesh->getPosY(i, j, k);
        FLOAT posZ = mesh->getPosZ(i, j, k);

        file << posX << " " << posY << " " << posZ << "\n";
      }
    }
  }

  file << "CELL_DATA " << cells << "\n";

  file << "SCALARS pressure float 1\n";
  file << "LOOKUP_TABLE default\n";

  for (auto &p : pressures) {
    file << p << "\n";
  }

  file << "VECTORS velocity float\n";

  for (int i = 0; i < cells; i++) {
    file << velocitiesX[i] << " "
         << velocitiesY[i] << " "
         << velocitiesZ[i] << "\n";
  }

  file.close();
}
