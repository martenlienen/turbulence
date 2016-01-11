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
  MPIIterator(FF& flowField, const Parameters& parameters, 
              std::vector<std::string> vec2D, std::vector<std::string> vec3D,
                std::function<void(FF& flowField, int, int, int, T&,std::vector<int>&)> apply)
      : Iterator<FF>(flowField, parameters),
        _flowField(flowField),
        _vec(parameters.geometry.dim==2?vec2D:vec3D),
        _veclocal(_vec.size(),-1),
        _size(0),
        _sizetotal(parameters.geometry.dim==2?
          flowField.getCellsX() *flowField.getCellsY()*_size:
          flowField.getCellsX() *flowField.getCellsY()*flowField.getCellsZ()*_size),
        _p(parameters),
        _apply(apply),
        _fname(parameters.vtk.prefix),
        _data(_sizetotal) {}

  void iterate();
  
  void fillTable();

 public:
  FF& _flowField;
  std::vector<std::string> _vec;
  std::vector<int> _veclocal;
  int _size;
  int _sizetotal;
  const Parameters& _p;
  std::function<void(FF& flowField, int, int, int, T&,std::vector<int>&)> _apply;
  std::string _fname;
  int counter = 0;
  double _scenario = 0;

 public:
  std::vector<T> _data;
  std::vector<std::vector<std::string>> _table2D { 
    {"p", "u", "v"},
    {"p", "u", "v", "nut"},
    {"p", "u", "v", "k", "epsilon","f1","f2","fmu","d","e"},
  };
  std::vector<std::vector<std::string>> _table3D { 
    {"p", "u", "v", "w"},
    {"p", "u", "v", "w", "nut"},
    {"p", "u", "v", "w", "k", "epsilon","f1","f2","fmu","d","e"},
  };
};

template <typename FF, typename T>
void MPIIterator<FF, T>::iterate() {

  
  if (this->_p.geometry.dim == 2) {
    for (int i = 0; i < _flowField.getCellsX(); i++) {
      for (int j = 0; j < _flowField.getCellsY(); j++) {
        
        this->_apply(this->_flowField, i, j, 0, this->_data[counter],this->_veclocal);
        counter += _size;
        
      }
    }
  } else {
    for (int i = 0; i <= _p.geometry.sizeX; i++) {
      for (int j = 0; j <= _p.geometry.sizeY; j++) {
        for (int k = 0; k <= _p.geometry.sizeZ; k++) {

          this->_apply(this->_flowField, i, j, k, this->_data[counter],this->_veclocal);
          counter += _size;
          
        }
      }
    }
  }
}

template <typename FF, typename T>
void MPIIterator<FF, T>::fillTable() {
  std::vector<std::vector<std::string>> table = 
      this->_p.geometry.dim == 2? this->_table2D : this->_table3D;
  
  std::vector<std::string> t = table[this->_scenario];
  
  for(unsigned long int j = 0; j < this->_vec.size();j++){
    for(unsigned long int k = 0; k < t.size();k++){
      if(this->_vec[j] == t[k]){
        this->_veclocal[j] = k;
      }
    }
  }

  this->_size = t.size();
  
  this->_sizetotal = this->_parameters.geometry.dim==2?
          this->_flowField.getCellsX() *this->_flowField.getCellsY()*_size:
          this->_flowField.getCellsX() *this->_flowField.getCellsY()*this->_flowField.getCellsZ()*_size;
  
  _data.resize(this->_sizetotal);
}

}

#endif  // RANGE_ITERATOR_H_
