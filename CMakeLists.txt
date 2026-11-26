cmake_minimum_required(VERSION 3.14)
project(NucMutInt LANGUAGES CXX)

set(CMAKE_CXX_STANDARD 17)
set(CMAKE_CXX_STANDARD_REQUIRED ON)
set(CMAKE_CXX_EXTENSIONS OFF)

find_package(Eigen3 3.3 REQUIRED NO_MODULE)

add_executable(NucMutInt
    src/HelixFromTmpl.cpp
    third_party/tinyxml2/tinyxml2.cpp
)

target_include_directories(NucMutInt PRIVATE
    ${CMAKE_CURRENT_SOURCE_DIR}/src
    ${CMAKE_CURRENT_SOURCE_DIR}/third_party/tinyxml2
)

target_link_libraries(NucMutInt PRIVATE Eigen3::Eigen)
