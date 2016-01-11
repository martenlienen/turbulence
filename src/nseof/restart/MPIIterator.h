#ifndef MPI_ITERATOR_H_
#define MPI_ITERATOR_H_

#include <functional>
#include <vector>

#include "../Iterators.h"
#include "../Point.h"

namespace nseof {

template <typename FF, typename T>
class MPIIterator : public Iterator<FF> {
 public:
  MPIIterator(FF& flowField, const Parameters& parameters, int size2, int size3,
                std::function<void(FF& flowField, int, int, int, T&)> apply)
      : Iterator<FF>(flowField, parameters),
        _flowField(flowField),
        _size(parameters.geometry.dim==2?size2:size3),
        _sizetotal(parameters.geometry.dim==2?
          flowField.getCellsX() *flowField.getCellsY()*_size:
          flowField.getCellsX() *flowField.getCellsY()*flowField.getCellsZ()*_size),
        _p(parameters),
        _apply(apply),
        _fname(parameters.vtk.prefix),
        _data(_sizetotal) {}

  void iterate();

 public:
  FF& _flowField;
  int _size;
  int _sizetotal;
  const Parameters& _p;
  std::function<void(FF& flowField, int, int, int, T&)> _apply;
  std::string _fname;

 public:
  std::vector<T> _data;
};

template <typename FF, typename T>
void MPIIterator<FF, T>::iterate() {

  int counter = 0;
  std::cout<<_data.size()<<std::endl;
  if (this->_p.geometry.dim == 2) {
    for (int i = 0; i < _flowField.getCellsX(); i++) {
      for (int j = 0; j < _flowField.getCellsY(); j++) {
        
        this->_apply(this->_flowField, i, j, 0, this->_data[counter]);
        counter += _size;
        
      }
    }
  } else {
    for (int i = 0; i <= _p.geometry.sizeX; i++) {
      for (int j = 0; j <= _p.geometry.sizeY; j++) {
        for (int k = 0; k <= _p.geometry.sizeZ; k++) {

          this->_apply(this->_flowField, i, j, k, this->_data[counter]);
          counter += _size;
          
        }
      }
    }
  }
}

}

#endif  // RANGE_ITERATOR_H_
