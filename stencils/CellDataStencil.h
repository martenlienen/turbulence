#ifndef _STENCILS_CELL_DATA_STENCIL_H
#define _STENCILS_CELL_DATA_STENCIL_H

#include <functional>
#include <memory>
#include <string>
#include <vector>

#include <vtk/SimpleCellData.h>

#include "../Stencil.h"
#include "../FlowField.h"
#include "../Parameters.h"

/**
 * Gather vtk::CellData from a flow field
 */
template <typename T>
class CellDataStencil : public FieldStencil<FlowField> {
 public:
  CellDataStencil(Parameters const& parameters, std::string dataName,
                  std::function<T(FlowField&, int, int)> apply2d,
                  std::function<T(FlowField&, int, int, int)> apply3d);

  void apply(FlowField& flowField, int i, int j);
  void apply(FlowField& flowField, int i, int j, int k);

  std::unique_ptr<vtk::SimpleCellData<T>> get();

  void reset();

 private:
  std::string dataName;
  std::vector<T> data;

  std::function<T(FlowField&, int, int)> apply2d;
  std::function<T(FlowField&, int, int, int)> apply3d;
};

// Include the implementation here because this is a template type

template <typename T>
CellDataStencil<T>::CellDataStencil(
    Parameters const& parameters, std::string dataName,
    std::function<T(FlowField&, int, int)> apply2d,
    std::function<T(FlowField&, int, int, int)> apply3d)
    : FieldStencil(parameters),
      dataName(dataName),
      data(),
      apply2d(apply2d),
      apply3d(apply3d) {}

template <typename T>
void CellDataStencil<T>::apply(FlowField& flowField, int i, int j) {
  this->data.push_back(this->apply2d(flowField, i, j));
}

template <typename T>
void CellDataStencil<T>::apply(FlowField& flowField, int i, int j, int k) {
  this->data.push_back(this->apply3d(flowField, i, j, k));
}

template <typename T>
std::unique_ptr<vtk::SimpleCellData<T>> CellDataStencil<T>::get() {
  return std::unique_ptr<vtk::SimpleCellData<T>>(
      new vtk::SimpleCellData<T>(this->dataName, this->data));
}

template <typename T>
void CellDataStencil<T>::reset() {
  this->data.clear();
}

#endif  // _STENCILS_CELL_DATA_STENCIL_H
