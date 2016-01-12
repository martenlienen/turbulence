#ifndef MPI_ITERATOR_READ_H_
#define MPI_ITERATOR_READ_H_

#include <functional>
#include <vector>
#include "MPIIterator.h"

#include "../Iterators.h"
#include "../Point.h"

namespace nseof {

template <typename FF, typename T>
class MPIIteratorRead : public MPIIterator<FF,T> {
 public:
  MPIIteratorRead(FF& flowField, const Parameters& parameters, 
                std::vector<std::string> vec2D, std::vector<std::string> vec3D,
                std::function<void(FF& flowField, int, int, int, T&,std::vector<int>&)> apply)
      : MPIIterator<FF,T>(flowField, parameters,vec2D, vec3D, apply) {}

  void iterate();

};

template <typename FF, typename T>
void MPIIteratorRead<FF, T>::iterate() {
  
  // is file existent
  if(!this->is_file_exist(this->_fname.c_str())){
    return;
  }
  
  // yes
  MPI_File fh, fh2;
  MPI_Status status, status2;
  int file_open_error, file_open_error2;
  MPI_Offset my_offset = 0;
  
  
  MPI_Barrier(PETSC_COMM_WORLD);
  
  
  // load simulation type
  if(this->_p.parallel.rank != 0){
    
    file_open_error2 = MPI_File_open(MPI_COMM_WORLD, this->_fname.c_str(), 
                                    MPI_MODE_RDONLY, MPI_INFO_NULL, &fh2);
    MPI_File_close(&fh2);
    

  } else{
    // try to open file
    file_open_error2 = MPI_File_open(MPI_COMM_WORLD, this->_fname.c_str(), 
                                    MPI_MODE_RDONLY, MPI_INFO_NULL, &fh2);
    
    if (file_open_error2 != MPI_SUCCESS) {
      // no success: nothing to load
    }else{
      // success: load data
      MPI_File_seek(fh2, 0, MPI_SEEK_SET);
      MPI_File_read(fh2, &(this->_scenario), 1, MPI_DOUBLE, &status2);
      MPI_File_close(&fh2);
      
    }
    
  }
  
  // broadcast simulation type
  MPI_Bcast(&(this->_scenario),1,MPI_DOUBLE,0, PETSC_COMM_WORLD);
  
  // fill in coversion table
  this->fillTable();
  
  // set offsets
  if(this->_p.parallel.rank != 0){
    my_offset = sizeof(T) * (this->_data.size() * 
      this->_p.parallel.rank + 1);
  } else {
    this->counter++;
  } 

  // read data
  // try to open file
  file_open_error = MPI_File_open(MPI_COMM_WORLD, this->_fname.c_str(), 
		                  MPI_MODE_RDONLY, MPI_INFO_NULL, &fh);
  
  if (file_open_error != MPI_SUCCESS) {
    // no success: nothing to load
  }else{
    // success: load data
    MPI_File_seek(fh, my_offset, MPI_SEEK_SET);
    MPI_File_read(fh, this->_data.data(), this->_sizetotal, MPI_DOUBLE, &status);
    MPI_File_close(&fh);
    
    // write data into flowfield
    MPIIterator<FF,T>::iterate();
    
  }
}

}

#endif  // RANGE_ITERATOR_H_
