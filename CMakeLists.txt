cmake_minimum_required(VERSION 3.11)

# add Eigen as submodule
# git submodule add https://gitlab.com/libeigen/eigen

project(EqSolver)

add_subdirectory(googletest)



add_subdirectory(src)

enable_testing()
add_subdirectory(test)




