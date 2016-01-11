#ifndef MPI_ITERATOR_WRITE_H_
#define MPI_ITERATOR_WRITE_H_

#include <functional>
#include <vector>
#include "MPIIterator.h"

#include "../Iterators.h"
#include "../Point.h"

namespace nseof {

template <typename FF, typename T>
class MPIIteratorWrite : public MPIIterator<FF,T> {
 public:
  MPIIteratorWrite(FF& flowField, const Parameters& parameters, 
                std::vector<std::string> vec2D, std::vector<std::string> vec3D,
                std::function<void(FF& flowField, int, int, int, T&,std::vector<int>&)> apply)
      : MPIIterator<FF,T>(flowField, parameters,vec2D, vec3D, apply) {
      
        std::string type = this->_p.simulation.type;
        
        if (type=="dns"){
          this->_scenario = 0;
        }else if(type=="aturb"){
          this->_scenario = 1;
        }else if(type=="keturb"){
          this->_scenario = 2;
        }
        
        this->fillTable();
      }

  void iterate();

};

template <typename FF, typename T>
void MPIIteratorWrite<FF, T>::iterate() {
  
  MPI_Offset my_offset = 0;
  
  if(this->_p.parallel.rank != 0){
    my_offset = sizeof(T) * (this->_data.size() * 
      this->_p.parallel.rank + 1);
  } else{
    // define type
    this->counter++;
    this->_data[0] = this->_scenario;
  }
  
  MPI_File fh;
  MPI_Status status;
  
  // load data from flowfield
  MPIIterator<FF,T>::iterate();
  
  // Wrire data to file
  MPI_File_open(MPI_COMM_WORLD, this->_fname.c_str(), 
		MPI_MODE_CREATE | MPI_MODE_RDWR, MPI_INFO_NULL, &fh);
  MPI_File_seek(fh, my_offset, MPI_SEEK_SET);
  MPI_File_write(fh, this->_data.data(), this->_sizetotal, MPI_DOUBLE, &status);
  MPI_File_close(&fh);
  
  
}

}

#endif  // RANGE_ITERATOR_H_
