cmake_minimum_required(VERSION 3.6)
project(L_shape_weak)

set(CMAKE_CXX_STANDARD 11)
set(CMAKE_CXX_FLAGS "-O2")
set(SOURCE_FILES main.cpp spline.cpp spline.h)
add_executable(L_shape_weak ${SOURCE_FILES})