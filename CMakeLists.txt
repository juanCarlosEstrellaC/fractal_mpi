cmake_minimum_required(VERSION 3.29)
project(ejemplo02_fractal2)

set(CMAKE_CXX_STANDARD 17)

find_package(fmt)
find_package(SFML)
find_package(glfw3)


find_package(MPI)
if (MPI_FOUND)
    message(STATUS "MPI FOUD")
    message(STATUS "MPI Include path: ${MPI_CXX_VERSION}")
    message(STATUS "MPI Include path: ${MPI_INCLUDE_PATH}")
    message(STATUS "MPI Libraries: ${MPI_LIBRARIES}")
    message(STATUS "MPI Include path: ${MPI_CXX_COMPILER}")
    message(STATUS "MPI Include path: ${MPI_lINK_FLAGS}")
else ()
    message(FATAL_ERROR "MPI Not FOUND")
endif ()
include_directories(${MPI_INCLUDE_PATH})


add_executable(ejemplo02_fractal2 main.cpp
        fps_counter.cpp
        fps_counter.h)
target_link_libraries(ejemplo02_fractal2 fmt::fmt)
target_link_libraries(ejemplo02_fractal2 sfml::sfml)
target_link_libraries(ejemplo02_fractal2 glfw)
target_link_libraries(ejemplo02_fractal2 ${MPI_LIBRARIES})
