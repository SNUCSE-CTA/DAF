cmake_minimum_required(VERSION 3.0)
project(DAF)

find_package(Boost 1.71.0 REQUIRED)
include_directories(${Boost_INCLUDE_DIRS})

include_directories(${PROJECT_SOURCE_DIR})

file(GLOB SOURCES src/*.cc)

add_subdirectory(main)
