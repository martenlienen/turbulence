#ifndef _STENCILS_CELL_DATA_STENCIL_H
#define _STENCILS_CELL_DATA_STENCIL_H

#include <functional>
#include <memory>
#include <string>
#include <vector>

#include <vtk/SimpleCellData.h>

#include "../Stencil.h"
#include "../Parameters.h"

namespace nseof {

/**
 * Gather vtk::CellData from a flow field
 */
template <typename T, typename FF>
class CellDataStencil : public FieldStencil<FF> {
 public:
  CellDataStencil(Parameters const& parameters, std::string dataName,
                  std::function<T(FF&, int, int)> apply2d,
                  std::function<T(FF&, int, int, int)> apply3d);

  void apply(FF& flowField, int i, int j);
  void apply(FF& flowField, int i, int j, int k);

  std::unique_ptr<vtk::SimpleCellData<T>> get();

  void reset();

 private:
  std::string dataName;
  std::vector<T> data;

  std::function<T(FF&, int, int)> apply2d;
  std::function<T(FF&, int, int, int)> apply3d;
};

// Include the implementation here because this is a template type

template <typename T, typename FF>
CellDataStencil<T, FF>::CellDataStencil(
    Parameters const& parameters, std::string dataName,
    std::function<T(FF&, int, int)> apply2d,
    std::function<T(FF&, int, int, int)> apply3d)
    : FieldStencil<FF>(parameters),
      dataName(dataName),
      data(),
      apply2d(apply2d),
      apply3d(apply3d) {}

template <typename T, typename FF>
void CellDataStencil<T, FF>::apply(FF& flowField, int i, int j) {
  if (flowField.getFlags().getValue(i, j) & OBSTACLE_SELF) {
    this->data.push_back(0);
  } else {
    this->data.push_back(this->apply2d(flowField, i, j));
  }
}

template <typename T, typename FF>
void CellDataStencil<T, FF>::apply(FF& flowField, int i, int j, int k) {
  if (flowField.getFlags().getValue(i, j, k) & OBSTACLE_SELF) {
    this->data.push_back(0);
  } else {
    this->data.push_back(this->apply3d(flowField, i, j, k));
  }
}

template <typename T, typename FF>
std::unique_ptr<vtk::SimpleCellData<T>> CellDataStencil<T, FF>::get() {
  return std::unique_ptr<vtk::SimpleCellData<T>>(
      new vtk::SimpleCellData<T>(this->dataName, this->data));
}

template <typename T, typename FF>
void CellDataStencil<T, FF>::reset() {
  this->data.clear();
}
}

#endif  // _STENCILS_CELL_DATA_STENCIL_H
