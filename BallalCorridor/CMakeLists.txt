
cmake_minimum_required(VERSION 3.10)

project(BallalCorridor LANGUAGES CXX)

set(CMAKE_CXX_STANDARD 20)

find_package(drake CONFIG REQUIRED PATHS /opt/drake)



add_library(BallalCorridor
    Ball.h
    Ball.cc
    GenerateBall.h
    BallCorridor.h
    GenerateBallCorridor.h
    LagrangePolynomial.h
    )

target_link_libraries(BallalCorridor drake::drake)

set_property(TARGET BallalCorridor PROPERTY CXX_STANDARD 20)

target_include_directories(BallalCorridor PUBLIC
  "${CMAKE_CURRENT_SOURCE_DIR}",
  "${CMAKE_CURRENT_SOURCE_DIR}/BallalCorridor"
)
