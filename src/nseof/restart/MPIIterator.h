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
  MPIIterator(FF& flowField, const Parameters& parameters, std::string fname,
              std::vector<std::string> vec2D, std::vector<std::string> vec3D,
              std::function<void(FF& flowField, int, int, int, T&,
                                 std::vector<int>&)> apply2D,
              std::function<void(FF& flowField, int, int, int, T&,
                                 std::vector<int>&)> apply3D)
      : Iterator<FF>(flowField, parameters),
        _flowField(flowField),
        _infocells(1),
        _vec(parameters.geometry.dim == 2 ? vec2D : vec3D),
        _veclocal(_vec.size(), -1),
        _size(0),
        _sizetotal(0),
        _data(0),
        _p(parameters),
        _apply(parameters.geometry.dim == 2 ? apply2D : apply3D),
        _fname(fname) {}

  void iterate();

  void fillTable();

 public:
  FF& _flowField;
  // the first cells in the vector of rank 0 are reserved for informations of
  // the simulation, at the moment:
  // [0] simulation type
  // at a later time possible:
  // [1] final time
  // [2] dimension
  // [3] ...
  int _infocells;
  std::vector<std::string> _vec;
  std::vector<int> _veclocal;
  int _size;
  int _sizetotal;
  std::vector<T> _data;
  const Parameters& _p;
  std::function<void(FF& flowField, int, int, int, T&, std::vector<int>&)>
      _apply;
  std::string _fname;
  int counter = 0;
  double _scenario = 0;

  std::vector<std::vector<std::string>> _table2D{
      {"p", "u", "v"},
      {"p", "u", "v", "nut"},
      {"p", "u", "v", "k", "epsilon", "f1", "f2", "fmu", "d", "e"},
  };
  std::vector<std::vector<std::string>> _table3D{
      {"p", "u", "v", "w"},
      {"p", "u", "v", "w", "nut"},
      {"p", "u", "v", "w", "k", "epsilon", "f1", "f2", "fmu", "d", "e"},
  };

  bool is_file_exist(const char* fileName) {
    std::ifstream infile(fileName);
    return infile.good();
  }
};

template <typename FF, typename T>
void MPIIterator<FF, T>::iterate() {
  if (this->_p.geometry.dim == 2) {
    for (int i = 0; i < _flowField.getCellsX(); i++) {
      for (int j = 0; j < _flowField.getCellsY(); j++) {
        this->_apply(this->_flowField, i, j, 0, this->_data[counter],
                     this->_veclocal);
        counter += _size;
      }
    }
  } else {
    for (int i = 0; i < _flowField.getCellsX(); i++) {
      for (int j = 0; j < _flowField.getCellsY(); j++) {
        for (int k = 0; k < _flowField.getCellsZ(); k++) {
          this->_apply(this->_flowField, i, j, k, this->_data[counter],
                       this->_veclocal);
          counter += _size;
        }
      }
    }
  }
}

template <typename FF, typename T>
void MPIIterator<FF, T>::fillTable() {
  // create conversion table:
  // format as defined in simulation -> format as in file
  std::vector<std::vector<std::string>> table =
      this->_p.geometry.dim == 2 ? this->_table2D : this->_table3D;

  std::vector<std::string> t = table[this->_scenario];

  for (unsigned long int j = 0; j < this->_vec.size(); j++) {
    for (unsigned long int k = 0; k < t.size(); k++) {
      if (this->_vec[j] == t[k]) {
        this->_veclocal[j] = k;
      }
    }
  }

  // parameters to save for each cell in FF
  this->_size = t.size();

  // new size of vector
  this->_sizetotal =
      this->_parameters.geometry.dim == 2
          ? this->_flowField.getCellsX() * this->_flowField.getCellsY() * _size
          : this->_flowField.getCellsX() * this->_flowField.getCellsY() *
                this->_flowField.getCellsZ() * _size;

  if (this->_p.parallel.rank == 0) {
    // in the first cell of binary: save type
    _sizetotal += _infocells;
  }

  // set size of vector
  _data.resize(this->_sizetotal);
}
}

#endif  // RANGE_ITERATOR_H_
