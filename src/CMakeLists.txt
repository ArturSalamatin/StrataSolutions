cmake_minimum_required(VERSION 3.11)


set(project "main")

add_executable(${project} main.cpp)

target_include_directories(
    ${project} PRIVATE 
      "${PROJECT_SOURCE_DIR}/../eigen"
  )