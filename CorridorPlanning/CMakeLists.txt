
cmake_minimum_required(VERSION 3.10)

project(CorridorPlanning LANGUAGES CXX)

set(CMAKE_CXX_STANDARD 20)

find_package(drake CONFIG REQUIRED PATHS /opt/drake)



add_library(CorridorPlanning
    CorridorPlanningProblem.h
    CorridorPlanningProblem.cc
    TimeVaryingProjrctionToBall.h
)

target_link_libraries(CorridorPlanning drake::drake)

set_property(TARGET CorridorPlanning PROPERTY CXX_STANDARD 20)

target_include_directories(CorridorPlanning PUBLIC
  "${CMAKE_CURRENT_SOURCE_DIR}",
  "${CMAKE_CURRENT_SOURCE_DIR}/BallalCorridor"
  "${CMAKE_CURRENT_SOURCE_DIR}/CorridorPlanning"
)
