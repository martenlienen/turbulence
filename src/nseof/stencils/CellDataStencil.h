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

template <typename T>
struct NullValue {};

template <>
struct NullValue<float> {
  static const float value;
};

const float NullValue<float>::value = 0.0;

template <>
struct NullValue<double> {
  static const double value;
};

const double NullValue<double>::value = 0.0;

template <>
struct NullValue<std::vector<double>> {
  static const std::vector<double> value;
};

const std::vector<double> NullValue<std::vector<double>>::value = {0.0, 0.0,
                                                                   0.0};

template <>
struct NullValue<std::vector<float>> {
  static const std::vector<float> value;
};

const std::vector<float> NullValue<std::vector<float>>::value = {0.0, 0.0,
                                                                   0.0};


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
    this->data.push_back(NullValue<T>::value);
  } else {
    this->data.push_back(this->apply2d(flowField, i, j));
  }
}

template <typename T, typename FF>
void CellDataStencil<T, FF>::apply(FF& flowField, int i, int j, int k) {
  if (flowField.getFlags().getValue(i, j, k) & OBSTACLE_SELF) {
    this->data.push_back(NullValue<T>::value);
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
