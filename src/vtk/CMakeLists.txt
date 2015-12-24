cmake_minimum_required(VERSION 2.8)
project("vtk")

add_compile_options(-std=c++14)

set(CMAKE_CXX_FLAGS_DEBUG
  "${CMAKE_CXX_FLAGS_DEBUG} -Wall -Werror -Wpedantic -Wextra")

set(SRCS File.cpp Dataset.cpp CellData.cpp)

add_library(vtk ${SRCS})
