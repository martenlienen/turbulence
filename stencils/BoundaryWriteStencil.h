#include <functional>
#include <utility>
#include <vector>

#include "../Stencil.h"

template <typename T, class FF>
class BoundaryWriteStencil : public BoundaryStencil<FF> {
 public:
  std::vector<T> leftData;
  std::vector<T> rightData;
  std::vector<T> topData;
  std::vector<T> bottomData;
  std::vector<T> frontData;
  std::vector<T> backData;

  BoundaryWriteStencil(const Parameters& parameters,
                       std::function<void(FF&, int, int, T)> apply2d,
                       std::function<void(FF&, int, int, int, T)> apply3d)
      : BoundaryStencil<FF>(parameters), apply2d(apply2d), apply3d(apply3d) {}
  ~BoundaryWriteStencil() {}

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
  std::function<void(FF&, int, int, T)> apply2d;
  std::function<void(FF&, int, int, int, T)> apply3d;
};

template <typename T, typename FF>
void BoundaryWriteStencil<T, FF>::applyLeftWall(FF& flowField, int i, int j) {
  this->apply2d(flowField, i, j, this->leftData.back());
  this->leftData.pop_back();
}

template <typename T, typename FF>
void BoundaryWriteStencil<T, FF>::applyRightWall(FF& flowField, int i, int j) {
  this->apply2d(flowField, i, j, this->rightData.back());
  this->rightData.pop_back();
}

template <typename T, typename FF>
void BoundaryWriteStencil<T, FF>::applyTopWall(FF& flowField, int i, int j) {
  this->apply2d(flowField, i, j, this->topData.back());
  this->topData.pop_back();
}

template <typename T, typename FF>
void BoundaryWriteStencil<T, FF>::applyBottomWall(FF& flowField, int i, int j) {
  this->apply2d(flowField, i, j, this->bottomData.back());
  this->bottomData.pop_back();
}

template <typename T, typename FF>
void BoundaryWriteStencil<T, FF>::applyLeftWall(FF& flowField, int i, int j,
                                                int k) {
  this->apply3d(flowField, i, j, k, this->leftData.back());
  this->leftData.pop_back();
}

template <typename T, typename FF>
void BoundaryWriteStencil<T, FF>::applyRightWall(FF& flowField, int i, int j,
                                                 int k) {
  this->apply3d(flowField, i, j, k, this->rightData.back());
  this->rightData.pop_back();
}

template <typename T, typename FF>
void BoundaryWriteStencil<T, FF>::applyTopWall(FF& flowField, int i, int j,
                                               int k) {
  this->apply3d(flowField, i, j, k, this->topData.back());
  this->topData.pop_back();
}

template <typename T, typename FF>
void BoundaryWriteStencil<T, FF>::applyBottomWall(FF& flowField, int i, int j,
                                                  int k) {
  this->apply3d(flowField, i, j, k, this->bottomData.back());
  this->bottomData.pop_back();
}

template <typename T, typename FF>
void BoundaryWriteStencil<T, FF>::applyFrontWall(FF& flowField, int i, int j,
                                                 int k) {
  this->apply3d(flowField, i, j, k, this->frontData.back());
  this->frontData.pop_back();
}

template <typename T, typename FF>
void BoundaryWriteStencil<T, FF>::applyBackWall(FF& flowField, int i, int j,
                                                int k) {
  this->apply3d(flowField, i, j, k, this->backData.back());
  this->backData.pop_back();
}
