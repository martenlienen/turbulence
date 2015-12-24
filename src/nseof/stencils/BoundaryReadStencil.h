#include <functional>
#include <vector>

#include "../Stencil.h"

namespace nseof {

template <typename T, class FF>
class BoundaryReadStencil : public BoundaryStencil<FF> {
 public:
  std::vector<T> leftData;
  std::vector<T> rightData;
  std::vector<T> topData;
  std::vector<T> bottomData;
  std::vector<T> frontData;
  std::vector<T> backData;

  BoundaryReadStencil(const Parameters& parameters,
                      std::function<T(FF&, int, int)> apply2d,
                      std::function<T(FF&, int, int, int)> apply3d)
      : BoundaryStencil<FF>(parameters), apply2d(apply2d), apply3d(apply3d) {}
  ~BoundaryReadStencil() {}

  void applyLeftWall(FF& flowField, int i, int j);

  void applyRightWall(FF& flowField, int i, int j);

  void applyBottomWall(FF& flowField, int i, int j);

  void applyTopWall(FF& flowField, int i, int j);

  void applyLeftWall(FF& flowField, int i, int j, int k);

  void applyRightWall(FF& flowField, int i, int j, int k);

  void applyTopWall(FF& flowField, int i, int j, int k);

  void applyBottomWall(FF& flowField, int i, int j, int k);

  void applyFrontWall(FF& flowField, int i, int j, int k);

  void applyBackWall(FF& flowField, int i, int j, int k);

 private:
  std::function<T(FF&, int, int)> apply2d;
  std::function<T(FF&, int, int, int)> apply3d;
};

// Include the implementation here because this is a template type

template <typename T, typename FF>
void BoundaryReadStencil<T, FF>::applyLeftWall(FF& flowField, int i, int j) {
  this->leftData.push_back(this->apply2d(flowField, i, j));
}

template <typename T, typename FF>
void BoundaryReadStencil<T, FF>::applyRightWall(FF& flowField, int i, int j) {
  this->rightData.push_back(this->apply2d(flowField, i, j));
}

template <typename T, typename FF>
void BoundaryReadStencil<T, FF>::applyTopWall(FF& flowField, int i, int j) {
  this->topData.push_back(this->apply2d(flowField, i, j));
}

template <typename T, typename FF>
void BoundaryReadStencil<T, FF>::applyBottomWall(FF& flowField, int i, int j) {
  this->bottomData.push_back(this->apply2d(flowField, i, j));
}

template <typename T, typename FF>
void BoundaryReadStencil<T, FF>::applyLeftWall(FF& flowField, int i, int j,
                                               int k) {
  this->leftData.push_back(this->apply3d(flowField, i, j, k));
}

template <typename T, typename FF>
void BoundaryReadStencil<T, FF>::applyRightWall(FF& flowField, int i, int j,
                                                int k) {
  this->rightData.push_back(this->apply3d(flowField, i, j, k));
}

template <typename T, typename FF>
void BoundaryReadStencil<T, FF>::applyTopWall(FF& flowField, int i, int j,
                                              int k) {
  this->topData.push_back(this->apply3d(flowField, i, j, k));
}

template <typename T, typename FF>
void BoundaryReadStencil<T, FF>::applyBottomWall(FF& flowField, int i, int j,
                                                 int k) {
  this->bottomData.push_back(this->apply3d(flowField, i, j, k));
}

template <typename T, typename FF>
void BoundaryReadStencil<T, FF>::applyFrontWall(FF& flowField, int i, int j,
                                                int k) {
  this->frontData.push_back(this->apply3d(flowField, i, j, k));
}

template <typename T, typename FF>
void BoundaryReadStencil<T, FF>::applyBackWall(FF& flowField, int i, int j,
                                               int k) {
  this->backData.push_back(this->apply3d(flowField, i, j, k));
}
}
