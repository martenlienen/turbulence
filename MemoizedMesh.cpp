#include "MemoizedMesh.h"

#include <utility>

/**
 * We add 2 everywhere because for example the BFStepInitStencil accesses cells
 * one step outside the local size.
 */
MemoizedMesh::MemoizedMesh(Parameters& params, std::unique_ptr<Meshsize> mesh)
    : MemoizedMesh(
          std::move(mesh), params.parallel.localSize[0] + 3 + 2,
          params.parallel.localSize[1] + 3 + 2,
          (params.geometry.dim == 3 ? params.parallel.localSize[2] : 0) + 3 + 2,
          params.parallel.localSize[0] + 4 + 2,
          params.parallel.localSize[1] + 4 + 2,
          (params.geometry.dim == 3 ? params.parallel.localSize[2] : 0) + 4 +
              2) {}

MemoizedMesh::MemoizedMesh(std::unique_ptr<Meshsize> mesh, int cellsX,
                           int cellsY, int cellsZ, int pointsX, int pointsY,
                           int pointsZ)
    : mesh(std::move(mesh)),
      cellsX(cellsX),
      cellsY(cellsY),
      cellsZ(cellsZ),
      pointsX(pointsX),
      pointsY(pointsY),
      pointsZ(pointsZ),
      dx(cellsX,
         std::vector<std::vector<FLOAT>>(cellsY, std::vector<FLOAT>(cellsZ))),
      dy(cellsX,
         std::vector<std::vector<FLOAT>>(cellsY, std::vector<FLOAT>(cellsZ))),
      dz(cellsX,
         std::vector<std::vector<FLOAT>>(cellsY, std::vector<FLOAT>(cellsZ))),
      posX(pointsX, std::vector<std::vector<FLOAT>>(
                        pointsY, std::vector<FLOAT>(pointsZ))),
      posY(pointsX, std::vector<std::vector<FLOAT>>(
                        pointsY, std::vector<FLOAT>(pointsZ))),
      posZ(pointsX, std::vector<std::vector<FLOAT>>(
                        pointsY, std::vector<FLOAT>(pointsZ))) {
  for (int i = 0; i < cellsX; i++) {
    for (int j = 0; j < cellsY; j++) {
      for (int k = 0; k < cellsZ; k++) {
        this->dx[i][j][k] = this->mesh->getDx(i - 1, j - 1, k - 1);
        this->dy[i][j][k] = this->mesh->getDy(i - 1, j - 1, k - 1);
        this->dz[i][j][k] = this->mesh->getDz(i - 1, j - 1, k - 1);
      }
    }
  }

  for (int i = 0; i < pointsX; i++) {
    for (int j = 0; j < pointsY; j++) {
      for (int k = 0; k < pointsZ; k++) {
        this->posX[i][j][k] = this->mesh->getPosX(i - 1, j - 1, k - 1);
        this->posY[i][j][k] = this->mesh->getPosY(i - 1, j - 1, k - 1);
        this->posZ[i][j][k] = this->mesh->getPosZ(i - 1, j - 1, k - 1);
      }
    }
  }
}
#include <iostream>
FLOAT MemoizedMesh::getDx(int i, int j) const {
  if (this->between(i, -1, this->cellsX - 2) &&
      this->between(j, -1, this->cellsY - 2)) {
    return this->dx[i + 1][j + 1][1];
  } else {
    return this->mesh->getDx(i, j, 0);
  }
}

FLOAT MemoizedMesh::getDy(int i, int j) const {
  if (this->between(i, -1, this->cellsX - 2) &&
      this->between(j, -1, this->cellsY - 2)) {
    return this->dy[i + 1][j + 1][1];
  } else {
    return this->mesh->getDx(i, j, 0);
  }
}

FLOAT MemoizedMesh::getDx(int i, int j, int k) const {
  if (this->between(i, -1, this->cellsX - 2) &&
      this->between(j, -1, this->cellsY - 2) &&
      this->between(k, -1, this->cellsZ)) {
    return this->dx[i + 1][j + 1][k + 1];
  } else {
    return this->mesh->getDx(i, j, 0);
  }
}

FLOAT MemoizedMesh::getDy(int i, int j, int k) const {
  if (this->between(i, -1, this->cellsX - 2) &&
      this->between(j, -1, this->cellsY - 2) &&
      this->between(k, -1, this->cellsZ)) {
    return this->dy[i + 1][j + 1][k + 1];
  } else {
    return this->mesh->getDx(i, j, 0);
  }
}

FLOAT MemoizedMesh::getDz(int i, int j, int k) const {
  if (this->between(i, -1, this->cellsX - 2) &&
      this->between(j, -1, this->cellsY - 2) &&
      this->between(k, -1, this->cellsZ)) {
    return this->dz[i + 1][j + 1][k + 1];
  } else {
    return this->mesh->getDx(i, j, 0);
  }
}

FLOAT MemoizedMesh::getPosX(int i, int j, int k) const {
  if (this->between(i, -1, this->pointsX - 2) &&
      this->between(j, -1, this->pointsY - 2) &&
      this->between(k, -1, this->pointsZ)) {
    return this->posX[i + 1][j + 1][k + 1];
  } else {
    return this->mesh->getDx(i, j, 0);
  }
}

FLOAT MemoizedMesh::getPosY(int i, int j, int k) const {
  if (this->between(i, -1, this->pointsX - 2) &&
      this->between(j, -1, this->pointsY - 2) &&
      this->between(k, -1, this->pointsZ)) {
    return this->posY[i + 1][j + 1][k + 1];
  } else {
    return this->mesh->getDx(i, j, 0);
  }
}

FLOAT MemoizedMesh::getPosZ(int i, int j, int k) const {
  if (this->between(i, -1, this->pointsX - 2) &&
      this->between(j, -1, this->pointsY - 2) &&
      this->between(k, -1, this->pointsZ)) {
    return this->posZ[i + 1][j + 1][k + 1];
  } else {
    return this->mesh->getDx(i, j, 0);
  }
}

FLOAT MemoizedMesh::getPosX(int i, int j) const {
  if (this->between(i, -1, this->pointsX - 2) &&
      this->between(j, -1, this->pointsY - 2)) {
    return this->posX[i + 1][j + 1][1];
  } else {
    return this->mesh->getDx(i, j, 0);
  }
}

FLOAT MemoizedMesh::getPosY(int i, int j) const {
  if (this->between(i, -1, this->pointsX - 2) &&
      this->between(j, -1, this->pointsY - 2)) {
    return this->posY[i + 1][j + 1][1];
  } else {
    return this->mesh->getDx(i, j, 0);
  }
}

FLOAT MemoizedMesh::getDxMin() const { return this->mesh->getDxMin(); }

FLOAT MemoizedMesh::getDyMin() const { return this->mesh->getDyMin(); }

FLOAT MemoizedMesh::getDzMin() const { return this->mesh->getDyMin(); }

bool MemoizedMesh::between(int x, int low, int high) const {
  return x >= low && x <= high;
}
