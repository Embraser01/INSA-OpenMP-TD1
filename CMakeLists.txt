cmake_minimum_required(VERSION 3.10)
project(td1)

set(CMAKE_CXX_STANDARD 14)

add_executable(td1 main.cpp)

message(STATUS "Checking OpenMP")
find_package(OpenMP)
IF (OPENMP_FOUND)
    message("Found OpenMP!")
    set(CMAKE_C_FLAGS " ${CMAKE_C_FLAGS} ${OpenMP_C_FLAGS}")
    set(CMAKE_CXX_FLAGS "${CMAKE_CXX_FLAGS} ${OpenMP_CXX_FLAGS}")
    set(CMAKE_SHARED_LINKER_FLAGS "${CMAKE_SHARED_LINKER_FLAGS} ${OpenMP_SHARED_LINKER_FLAGS}")
    set(CMAKE_EXE_LINKER_FLAGS "${CMAKE_EXE_LINKER_FLAGS} ${OpenMP_EXE_LINKER_FLAGS}")
ELSE ()
    message("Missed OpenMP! ")
ENDIF ()