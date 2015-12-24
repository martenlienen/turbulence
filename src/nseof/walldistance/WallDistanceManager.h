#ifndef _WALLDISTANCEMANAGER_H_
#define _WALLDISTANCEMANAGER_H_

#include <vector>

#include <ANN/ANN.h>

#include "../Parameters.h"

namespace nseof {

class WallDistanceManager {
 public:
  WallDistanceManager(const Parameters& parameters);

  virtual ~WallDistanceManager();

  void init();

  double query(int i, int j, int kk = 0);

 private:
  const Parameters& _parameters;
  std::vector<double> x;
  std::vector<double> y;
  std::vector<double> z;
  int _dim;
  std::string _scenario;
  int _sizex;
  int _sizey;
  int _sizez;

  // ANN fields
  // actual number of data points
  int nPts;

  // error bound
  double eps;

  // number of near neighbors
  int k;

  // data points
  ANNpointArray dataPts;

  // query point
  ANNpoint queryPt;

  // near neighbor indices
  ANNidxArray nnIdx;

  // near neighbor distances
  ANNdistArray dists;

  // search structure
  ANNkd_tree* kdTree;

  void initObstacles();
};
}

#endif  // _WALLDISTANCEMANAGER_H_
