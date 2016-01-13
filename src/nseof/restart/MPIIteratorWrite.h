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
  MPIIteratorWrite(FF& flowField, const Parameters& parameters, int size2, int size3,
                std::function<void(FF& flowField, int, int, int, T&)> apply)
      : MPIIterator<FF,T>(flowField, parameters,size2, size3, apply) {}

  void iterate();

};

template <typename FF, typename T>
void MPIIteratorWrite<FF, T>::iterate() {
  
  MPI_Offset my_offset = (sizeof(T) * this->_data.size()) * 
      this->_p.parallel.rank;
  MPI_File fh;
  MPI_Status status;
  
  // load data from flowfield
  MPIIterator<FF,T>::iterate();
  
  MPI_Barrier(MPI_COMM_WORLD);

  // Wrire data to file
  MPI_File_open(MPI_COMM_WORLD, this->_fname.c_str(), 
		MPI_MODE_CREATE | MPI_MODE_RDWR, MPI_INFO_NULL, &fh);
  MPI_File_seek(fh, my_offset, MPI_SEEK_SET);
  MPI_File_write(fh, this->_data.data(), this->_sizetotal, MPI_DOUBLE, &status);
  MPI_File_close(&fh);
  
  MPI_Barrier(MPI_COMM_WORLD);
  
}

}

#endif  // RANGE_ITERATOR_H_
