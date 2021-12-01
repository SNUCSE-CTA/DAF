cmake_minimum_required(VERSION 3.0)
project(DAF)

add_compile_options(-O3)
add_compile_options(-Wall)
add_compile_options(-std=c++17)
# set(CMAKE_CXX_FLAGS_DEBUG "${CMAKE_CXX_FLAGS_DEBUG} -fno-omit-frame-pointer -fsanitize=address")
# set(CMAKE_LINKER_FLAGS_DEBUG "${CMAKE_LINKER_FLAGS_DEBUG} -fno-omit-frame-pointer -fsanitize=address")

include_directories(${PROJECT_SOURCE_DIR})

file(GLOB SOURCES src/*.cc)

add_subdirectory(main)
