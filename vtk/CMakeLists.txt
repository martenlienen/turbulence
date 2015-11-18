cmake_minimum_required(VERSION 2.8)
project("vtk")

add_compile_options(-Wall -Werror -Wpedantic -Wextra -O3 -std=c++14)

set(SRCS File.cpp Dataset.cpp CellData.cpp)

add_library(vtk ${SRCS})
