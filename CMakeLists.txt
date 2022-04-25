cmake_minimum_required(VERSION 3.22)
project(TravelingSalesman C)

set(CMAKE_C_STANDARD 99)

add_executable(TravelingSalesman TravelingSalesmanMPI.c Util.c Util.h)
