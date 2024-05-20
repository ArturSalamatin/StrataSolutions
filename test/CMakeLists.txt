cmake_minimum_required(VERSION 3.11)


# add googletest as submodule
# <cd> git submodule add https://github.com/google/googletest  
# bind to a particular tag
# cd googletest
# git checkout v1.14.0


add_executable(grid_test grid_test.cpp)
target_include_directories(grid_test PRIVATE "${CMAKE_CURRENT_SOURCE_DIR}/../")
target_link_libraries(grid_test PRIVATE gtest gtest_main)

if(MSVC)
  # hack for windows compilation
  # or for compilation with MSVC
  # if MinGW -- may not be required...
  set(gtest_force_shared_crt ON CACHE BOOL "" FORCE)
endif()

add_test(grid_test grid_test)