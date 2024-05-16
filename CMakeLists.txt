cmake_minimum_required(VERSION 3.11)

# add Eigen as submodule
# git submodule add https://gitlab.com/libeigen/eigen

project(EqSolver)

# Sepcify standard ABOVE the call to add_executable
set(CMAKE_CXX_STANDARD 17)
set(CMAKE_CXX_EXTENSIONS OFF)
set(CMAKE_CXX_STANDARD_REQUIRED TRUE)

add_subdirectory(src)

enable_testing()
add_subdirectory(test)




