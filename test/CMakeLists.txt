cmake_minimum_required(VERSION 3.11)


# add googletest as submodule
# <cd> git submodule add https://github.com/google/googletest  
# bind to a particular tag
# cd googletest
# git checkout v1.14.0


macro(add_google_test name)
  set(TEST_EXECUTABLE ${name})
  set(TEST_NAME ${name})

  add_executable(
    ${TEST_NAME}  
      ${TEST_EXECUTABLE}.cpp
  )
  target_include_directories(
    ${TEST_NAME} PRIVATE 
      "${CMAKE_CURRENT_SOURCE_DIR}/../"
  )
  target_include_directories(
    ${TEST_NAME} PRIVATE 
      "${CMAKE_CURRENT_SOURCE_DIR}/../eigen"
  )
  target_link_libraries(
    ${TEST_NAME} PRIVATE
      gtest_main
  ) 

  # let the preprocessor know about the system name
  if(CMAKE_SYSTEM_NAME STREQUAL "Linux")
    target_compile_definitions(${TEST_NAME} PUBLIC "IS_LINUX")
  endif()
  if(CMAKE_SYSTEM_NAME STREQUAL "Darwin")
    target_compile_definitions(${TEST_NAME} PUBLIC "IS_MACOS")
  endif()
  if(CMAKE_SYSTEM_NAME STREQUAL "Windows")
    target_compile_definitions(${TEST_NAME} PUBLIC "IS_WINDOWS")
  endif()

  # compiler-specific options
  if(MSVC AND MSVC_VERSION GREATER 1400)
    target_compile_options(
      ${TEST_NAME} PRIVATE 
      /MP)
  else()
    target_compile_options(
      ${TEST_NAME} PRIVATE
        -Wall -Wextra -pedantic -Werror
    )
  endif()
  
  if(MSVC)
    # hack for windows compilation
    # or for compilation with MSVC
    # if MinGW -- may not be required...
    set(gtest_force_shared_crt ON CACHE BOOL "" FORCE)
  endif()

  # adding the test, finally
  add_test(
    NAME ${TEST_NAME} 
    COMMAND ${TEST_EXECUTABLE})

endmacro()




add_google_test(grid1D_test)
add_google_test(grid2D_test)
add_google_test(splitx_test)
add_google_test(splitting_solver_FactoryWithSource)
add_google_test(splitting_solver_SimpleProblem)
add_google_test(splitting_solver_SelfSimilarProblem)