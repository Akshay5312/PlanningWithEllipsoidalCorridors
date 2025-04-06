

cmake_minimum_required(VERSION 3.10)

project(CollisionObjects LANGUAGES CXX)

set(CMAKE_CXX_STANDARD 20)

option(DRAKE_PATH "The path for the drake installation" /opt/drake)

find_package(drake CONFIG REQUIRED PATHS /opt/drake)


# add_subdirectory(EllipsoidalCorridor)

# add_executable(CreateCorridors examples/CreateCorridors.cc)
# add_executable(CreateEllipsoid examples/CreateEllipsoid.cc)

# target_link_libraries(CreateCorridors PUBLIC drake::drake)
# target_link_libraries(CreateEllipsoid PUBLIC drake::drake)

add_subdirectory(BallalCorridor)

add_executable(GenerateLargeBall examples/GenerateLargeBall.cc)
target_link_libraries(GenerateLargeBall PUBLIC drake::drake)

add_executable(GenerateLargeBallCorridor examples/GenerateLargeBallCorridor.cc)
target_link_libraries(GenerateLargeBallCorridor PUBLIC drake::drake)