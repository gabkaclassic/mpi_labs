cmake_minimum_required(VERSION 3.26)
project(MPI_labs)

set(CMAKE_CXX_STANDARD 23)

find_package(MPI REQUIRED)

add_executable(MPI_labs main.cpp)

target_link_libraries(MPI_labs PUBLIC MPI::MPI_CXX)