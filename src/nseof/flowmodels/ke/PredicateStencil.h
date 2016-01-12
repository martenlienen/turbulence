#ifndef _NSEOF_FLOWMODELS_KE_PREDICATE_STENCIL_H_
#define _NSEOF_FLOWMODELS_KE_PREDICATE_STENCIL_H_

#include <string>
#include <fstream>
#include <sstream>
#include <iostream>
#include <functional>

#include "../../Definitions.h"
#include "../../Parameters.h"
#include "../../Stencil.h"

namespace nseof {

namespace flowmodels {

namespace ke {

using namespace std;

/** TODO WS1: Stencil for writting VTK files
 *
 * When iterated with, creates a VTK file.
 */
template <class FF, class SS>
class PredicateStencil : public FieldStencil<FF> {
 public:
  /** Constructor
   *
   * @param prefix String with the prefix of the name of the VTK files
   */
  PredicateStencil(
      const Parameters& parameters, double resetValue,
      std::function<bool(SS, SS)> predicate,
      std::function<SS(int, int, int, const Parameters&, FF&)> func);

  ~PredicateStencil();

  /** 2D operation for one position
   *
   * @param flowField State of the flow field
   * @param i Position in the x direction
   * @param j Position in the y direction
   */
  void apply(FF& flowField, int i, int j);

  /** 3D operation for one position
   *
   * @param flowField State of the flow field
   * @param i Position in the x direction
   * @param j Position in the y direction
   * @param k Position in the z direction
   */
  void apply(FF& flowField, int i, int j, int k);

  SS getValue();

  void reset();

 private:
  double _resetValue;
  double _tempValue;
  double _value;

  std::function<bool(SS, SS)> _predicate;
  std::function<SS(int, int, int, const Parameters&, FF&)> _func;

  const Parameters& _p;
};

template <class FF, class SS>
PredicateStencil<FF, SS>::PredicateStencil(
    const Parameters& parameters, double resetValue,
    std::function<bool(SS, SS)> predicate,
    std::function<SS(int, int, int, const Parameters&, FF&)> func)
    : FieldStencil<FF>(parameters),
      _resetValue(resetValue),
      _tempValue(resetValue),
      _value(resetValue),
      _predicate(predicate),
      _func(func),
      _p(parameters) {}

template <class FF, class SS>
PredicateStencil<FF, SS>::~PredicateStencil() {}

template <class FF, class SS>
void PredicateStencil<FF, SS>::apply(FF& flowField, int i, int j) {
  apply(flowField, i, j, 0);
}

template <class FF, class SS>
void PredicateStencil<FF, SS>::apply(FF& flowField, int i, int j, int k) {
  _tempValue = _func(i, j, k, _p, flowField);
  _value = _predicate(_value, _tempValue) ? _value : _tempValue;
}

template <class FF, class SS>
void PredicateStencil<FF, SS>::reset() {
  _value = _resetValue;
}

template <class FF, class SS>
SS PredicateStencil<FF, SS>::getValue() {
  return _value;
}
}
}
}

#endif
