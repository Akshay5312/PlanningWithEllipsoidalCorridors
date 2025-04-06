
cmake_minimum_required(VERSION 3.10)

project(EllipsoidalCorridor LANGUAGES CXX)

set(CMAKE_CXX_STANDARD 20)

find_package(drake CONFIG REQUIRED PATHS /opt/drake)



add_library(EllipsoidalCorridor
    EllipsoidalCorridor.h
    EllipsoidalCorridor.cc
    )

target_link_libraries(EllipsoidalCorridor drake::drake)

set_property(TARGET EllipsoidalCorridor PROPERTY CXX_STANDARD 20)

target_include_directories(EllipsoidalCorridor PUBLIC
  "${CMAKE_CURRENT_SOURCE_DIR}",
  "${CMAKE_CURRENT_SOURCE_DIR}/EllipsoidalCorridor"
)
