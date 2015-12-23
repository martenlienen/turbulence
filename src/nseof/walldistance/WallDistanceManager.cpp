#include "WallDistanceManager.h"

#include <iostream>

WallDistanceManager::WallDistanceManager(const Parameters& parameters)
    : _parameters(parameters),
      x(0),
      y(0),
      z(0),
      _dim(parameters.geometry.dim),
      _scenario(_parameters.simulation.scenario),
      _sizex(parameters.geometry.sizeX),
      _sizey(parameters.geometry.sizeY),
      _sizez(parameters.geometry.sizeZ),
      nPts(0),
      eps(0),
      k(1) {}

WallDistanceManager::~WallDistanceManager() {
  delete[] nnIdx;  // clean things up
  delete[] dists;
  delete kdTree;
  annClose();  // done with ANN
}

void WallDistanceManager::initObstacles() {
  // Pressure Channel - 2D
  if (_dim == 2) {
    // bottom-line
    for (int i = 0; i < _sizex + 2; i++) {
      x.push_back(i);
      y.push_back(1);
    }

    // top-line
    for (int i = 0; i < _sizex + 2; i++) {
      x.push_back(i);
      y.push_back(_sizey + 2);
    }
  }

  // Pressure Channel - 3D
  if (_dim == 3) {
    // bottom-surface
    for (int i = 0; i < _sizex + 2; i++) {
      for (int j = 0; j < _sizez + 2; j++) {
        x.push_back(i);
        y.push_back(1);
        z.push_back(j);
      }
    }

    // top-surface
    for (int i = 0; i < _sizex + 2; i++) {
      for (int j = 0; j < _sizez + 2; j++) {
        x.push_back(i);
        y.push_back(_sizey + 2);
        z.push_back(j);
      }
    }

    // front-surface
    for (int i = 0; i < _sizex + 2; i++) {
      for (int j = 0; j < _sizey + 2; j++) {
        x.push_back(i);
        y.push_back(j);
        z.push_back(1);
      }
    }

    // back-surface
    for (int i = 0; i < _sizex + 2; i++) {
      for (int j = 0; j < _sizey + 2; j++) {
        x.push_back(i);
        y.push_back(j);
        z.push_back(_sizez + 2);
      }
    }
  }

  // Channel - 2D
  if (_scenario == "channel") {
    double lx = _parameters.bfStep.xRatio * _parameters.geometry.lengthX;
    double ly = _parameters.bfStep.yRatio * _parameters.geometry.lengthY;

    // is there any step?
    if (lx > 0.0 && ly > 0) {
      // 2D
      if (_dim == 2) {
        // step-top
        for (int i = 0; i < _sizex + 2; i++) {
          if (_parameters.meshsize->getPosCellX(i, 0, 0) < lx) {
            x.push_back(i);
            y.push_back(-ly);

          } else {
            break;
          }
        }
        // step-right
        for (int i = 0; i < _sizey + 2; i++) {
          if (_parameters.meshsize->getPosCellY(0, i, 0) < ly) {
            x.push_back(-lx);
            y.push_back(i);

          } else {
            break;
          }
        }
      }

      // 3D
      if (_dim == 3) {
        // step-top
        for (int i = 0; i < _sizex + 2; i++) {
          if (_parameters.meshsize->getPosCellX(i, 0, 0) < lx) {
            for (int j = 0; j < _sizez + 2; j++) {
              x.push_back(i);
              y.push_back(-ly);
              z.push_back(j);
            }
          } else {
            break;
          }
        }
        // step-right
        for (int i = 0; i < _sizey + 2; i++) {
          if (_parameters.meshsize->getPosCellY(0, i, 0) < ly) {
            for (int j = 0; j < _sizez + 2; j++) {
              x.push_back(-lx);
              y.push_back(i);
              z.push_back(j);
            }
          } else {
            break;
          }
        }
      }
    }
  }

  // Cavity - 2D
  if (_scenario == "cavity" && _dim == 2) {
    // left-line
    for (int i = 0; i < _sizex + 2; i++) {
      x.push_back(1);
      y.push_back(i);
    }

    // right-line
    for (int i = 0; i < _sizex + 2; i++) {
      x.push_back(_sizex + 2);
      y.push_back(i);
    }
  }

  // Cavity - 3D
  if (_scenario == "cavity" && _dim == 3) {
    // left-surface
    for (int i = 0; i < _sizez + 2; i++) {
      for (int j = 0; j < _sizey + 2; j++) {
        x.push_back(1);
        y.push_back(j);
        z.push_back(i);
      }
    }

    // right-surface
    for (int i = 0; i < _sizez + 2; i++) {
      for (int j = 0; j < _sizey + 2; j++) {
        x.push_back(_sizex + 2);
        y.push_back(j);
        z.push_back(i);
      }
    }
  }
}

void WallDistanceManager::init() {
  // 1st STEP: load all data points (obstacles)
  if (_scenario == "pressure-channel" || _scenario == "channel" ||
      _scenario == "cavity") {
    initObstacles();
  }

  if (_scenario == "free") {
    // TODO
    // an arbitrary geometry could be loaded at this position
    // i.e. load obstacle cells from a file
  }

  // 2nd STEP: init ANN (size and dim of k-D-tree)
  int maxPts = x.size();
  nPts = 0;                             // read data points
  queryPt = annAllocPt(_dim);           // allocate query point
  dataPts = annAllocPts(maxPts, _dim);  // allocate data points
  nnIdx = new ANNidx[k];                // allocate near neigh indices
  dists = new ANNdist[k];               // allocate near neighbor dists

  // 3rd STEP: insert coordinates of data points in kd-tree
  if (_dim == 2) {
    for (unsigned long int i = 0; i < x.size(); i++) {
      double xx =
          x[i] < 0.0 ? -x[i] : _parameters.meshsize->getPosCellX(x[i], 0, 0);
      double yy =
          y[i] < 0.0 ? -y[i] : _parameters.meshsize->getPosCellY(0, y[i], 0);
      dataPts[nPts][0] = xx;
      dataPts[nPts][1] = yy;
      nPts++;
    }
  } else {
    for (unsigned long int i = 0; i < x.size(); i++) {
      double xx =
          x[i] < 0.0 ? -x[i] : _parameters.meshsize->getPosCellX(x[i], 0, 0);
      double yy =
          y[i] < 0.0 ? -y[i] : _parameters.meshsize->getPosCellY(0, y[i], 0);
      double zz =
          z[i] < 0.0 ? -z[i] : _parameters.meshsize->getPosCellZ(0, 0, z[i]);
      dataPts[nPts][0] = xx;
      dataPts[nPts][1] = yy;
      dataPts[nPts][2] = zz;
      nPts++;
    }
  }

  kdTree = new ANNkd_tree(  // build search structure
      dataPts,              // the data points
      nPts,                 // number of points
      _dim);                // dimension of space
}

double WallDistanceManager::query(int i, int j, int kk) {
  // 4st STEP: create query point
  queryPt[0] = _parameters.meshsize->getPosX(i, j, kk) +
               _parameters.meshsize->getDx(i, j, kk) / 2;
  queryPt[1] = _parameters.meshsize->getPosY(i, j, kk) +
               _parameters.meshsize->getDy(i, j, kk) / 2;

  if (_dim == 3) {
    queryPt[2] = _parameters.meshsize->getPosZ(i, j, kk) +
                 _parameters.meshsize->getDz(i, j, kk) / 2;
  }

  // 5st STEP: search shortest distance between query point and data points
  kdTree->annkSearch(  // search
      queryPt,         // query point
      k,               // number of near neighbors
      nnIdx,           // nearest neighbors (returned)
      dists,           // distance (returned)
      eps);            // error bound

  for (int i = 0; i < k; i++) {
    dists[i] = sqrt(dists[i]);  // unsquare distance
  }

  return dists[0];
}
