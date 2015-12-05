#ifndef _MPI_COMMUNICATOR_H_
#define _MPI_COMMUNICATOR_H_

#include <array>
#include <cstring>
#include <functional>
#include <memory>
#include <vector>

#include <mpi.h>

#include "../FlowField.h"
#include "../Parameters.h"

#include "../iterators/BottomTopBoundaryIterator.h"
#include "../iterators/FrontBackBoundaryIterator.h"
#include "../iterators/LeftRightBoundaryIterator.h"
#include "../iterators/ParallelBoundaryIterator.h"

template <typename T>
struct MPIType {
  typedef char type;
  static const MPI_Datatype mpiType;
};

template <typename T>
const MPI_Datatype MPIType<T>::mpiType = MPI_CHAR;

template <>
struct MPIType<float> {
  typedef float type;
  static const MPI_Datatype mpiType;
};

const MPI_Datatype MPIType<float>::mpiType = MPI_FLOAT;

template <>
struct MPIType<double> {
  typedef double type;
  static const MPI_Datatype mpiType;
};

const MPI_Datatype MPIType<double>::mpiType = MPI_DOUBLE;

template <>
struct MPIType<std::array<double, 3>> {
  typedef double type;
  static const MPI_Datatype mpiType;
};

const MPI_Datatype MPIType<std::array<double, 3>>::mpiType = MPI_DOUBLE;

template <typename T, typename FF>
class MPICommunicator {
 public:
  MPICommunicator(FF& flowField, Parameters& parameters,
                  std::function<void(FF&, int, int, int, T&)> read,
                  std::function<void(FF&, int, int, int, T&)> write,
                  int lbf = 1)
      : parameters(parameters),
        read(read),
        write(write),
        // clang-format off
        lrReadIterator(flowField, parameters,  {+2      , +2      , +2      }, {-1, -1, -1}, this->read , 1  , lbf),
        lrWriteIterator(flowField, parameters, {+2 - lbf, +2      , +2      }, {+0, -1, -1}, this->write, lbf, 1  ),
        btReadIterator(flowField, parameters,  {+2 - lbf, +2      , +2      }, {+0, -1, -1}, this->read , 1  , lbf),
        btWriteIterator(flowField, parameters, {+2 - lbf, +2 - lbf, +2      }, {+0, +0, -1}, this->write, lbf, 1  ),
        fbReadIterator(flowField, parameters,  {+2 - lbf, +2 - lbf, +2      }, {+0, +0, -1}, this->read , 1  , lbf),
        fbWriteIterator(flowField, parameters, {+2 - lbf, +2 - lbf, +2 - lbf}, {+0, +0, +0}, this->write, lbf, 1  ) {}
  // clang-format on

  void communicate(FF& flowField);

 private:
  Parameters& parameters;
  std::function<void(FF&, int, int, int, T&)> read;
  std::function<void(FF&, int, int, int, T&)> write;
  LeftRightBoundaryIterator<FF, T> lrReadIterator;
  LeftRightBoundaryIterator<FF, T> lrWriteIterator;
  BottomTopBoundaryIterator<FF, T> btReadIterator;
  BottomTopBoundaryIterator<FF, T> btWriteIterator;
  FrontBackBoundaryIterator<FF, T> fbReadIterator;
  FrontBackBoundaryIterator<FF, T> fbWriteIterator;

  void exchangeValues(ParallelBoundaryIterator<FF, T>& lrReadIterator,
                      ParallelBoundaryIterator<FF, T>& lrWriteIterator, int lnb,
                      int rnb);
};

// Implementation for template type

#include <mpi.h>

#include "Iterators.h"

template <typename T, typename FF>
void MPICommunicator<T, FF>::communicate(FF& flowField) {
  Parameters& parameters = this->parameters;

  this->exchangeValues(this->lrReadIterator, this->lrWriteIterator,
                       parameters.parallel.leftNb, parameters.parallel.rightNb);
  this->exchangeValues(this->btReadIterator, this->btWriteIterator,
                       parameters.parallel.bottomNb, parameters.parallel.topNb);
  this->exchangeValues(this->fbReadIterator, this->fbWriteIterator,
                       parameters.parallel.frontNb, parameters.parallel.backNb);
}

template <typename T, typename FF>
void MPICommunicator<T, FF>::exchangeValues(
    ParallelBoundaryIterator<FF, T>& readIterator,
    ParallelBoundaryIterator<FF, T>& writeIterator, int lnb, int rnb) {
  MPI_Datatype type = MPIType<T>::mpiType;
  MPI_Status status;

  // The raw data type that is transmitted by MPI (e.g. int or double)
  typedef typename MPIType<T>::type raw;

  std::vector<T>& lData = readIterator.first.data;
  std::vector<T>& rData = readIterator.second.data;

  // Number of primitive values (e.g. doubles) that we are sending left
  // respectively right
  int lNumSend = 0;
  int lNumRecv = 0;
  int rNumSend = 0;
  int rNumRecv = 0;

  if (lnb >= 0) {
    lNumSend = lData.size() * sizeof(T) / sizeof(raw);
    lNumRecv = rData.size() * sizeof(T) / sizeof(raw);

    readIterator.first.iterate();
  }

  if (rnb >= 0) {
    rNumSend = rData.size() * sizeof(T) / sizeof(raw);
    rNumRecv = lData.size() * sizeof(T) / sizeof(raw);

    readIterator.second.iterate();
  }

  // Buffers to hold the received values. Note that we receive as many bytes
  // from the left as we send to the right and vice versa.
  std::unique_ptr<raw> lBuf(new raw[lNumRecv]);
  std::unique_ptr<raw> rBuf(new raw[rNumRecv]);

  // clang-format off
  MPI_Sendrecv(lData.data(), lNumSend, type, lnb, 0,
               rBuf.get(), rNumRecv, type, rnb, 0,
               MPI_COMM_WORLD, &status);

  MPI_Sendrecv(rData.data(), rNumSend, type, rnb, 0,
               lBuf.get(), lNumRecv, type, lnb, 0,
               MPI_COMM_WORLD, &status);
  // clang-format on

  if (lnb >= 0) {
    // We use memcpy here because we really need to do a byte-by-byte copy of
    // the memory. std::copy would try to do actual type conversions which lead
    // to errors because there is no constructor to create for example an
    // std::array from a char.
    std::memcpy(writeIterator.first.data.data(), lBuf.get(),
                lNumRecv * sizeof(raw));

    writeIterator.first.iterate();
  }

  if (rnb >= 0) {
    std::memcpy(writeIterator.second.data.data(), rBuf.get(),
                rNumRecv * sizeof(raw));

    writeIterator.second.iterate();
  }
}

#endif  // _MPI_COMMUNICATOR_H_
