#ifndef _MEMOIZED_MESH_H_
#define _MEMOIZED_MESH_H_

#include <memory>
#include <vector>

#include "Meshsize.h"
#include "Parameters.h"

class MemoizedMesh : public Meshsize {
 public:
  MemoizedMesh(Parameters& params, std::unique_ptr<Meshsize> mesh);

  FLOAT getDx(int i, int j) const;
  FLOAT getDy(int i, int j) const;

  FLOAT getDx(int i, int j, int k) const;
  FLOAT getDy(int i, int j, int k) const;
  FLOAT getDz(int i, int j, int k) const;

  FLOAT getPosX(int i, int j, int k) const;
  FLOAT getPosY(int i, int j, int k) const;
  FLOAT getPosZ(int i, int j, int k) const;

  FLOAT getPosX(int i, int j) const;
  FLOAT getPosY(int i, int j) const;

  FLOAT getDxMin() const;
  FLOAT getDyMin() const;
  FLOAT getDzMin() const;

 private:
  std::unique_ptr<Meshsize> mesh;

  std::vector<std::vector<std::vector<FLOAT>>> dx;
  std::vector<std::vector<std::vector<FLOAT>>> dy;
  std::vector<std::vector<std::vector<FLOAT>>> dz;

  std::vector<std::vector<std::vector<FLOAT>>> posX;
  std::vector<std::vector<std::vector<FLOAT>>> posY;
  std::vector<std::vector<std::vector<FLOAT>>> posZ;

  MemoizedMesh(std::unique_ptr<Meshsize> mesh, int cellsX, int cellsY,
               int cellsZ, int pointsX, int pointsY, int pointsZ);
};

#endif  // _MEMOIZED_MESH_H_
