cmake_minimum_required(VERSION 2.8)
project("ANN" CXX)

include_directories(include)
include_directories(src)

file(GLOB SRCS src/*.cpp)

add_library(ANN ${SRCS})
