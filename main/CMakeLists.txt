add_compile_options(-O3)
add_compile_options(-Wall)
add_compile_options(-std=c++17)
add_compile_options(-DNDEBUG)
# set(CMAKE_CXX_FLAGS_DEBUG "${CMAKE_CXX_FLAGS_DEBUG} -fno-omit-frame-pointer -fsanitize=address")
# set(CMAKE_LINKER_FLAGS_DEBUG "${CMAKE_LINKER_FLAGS_DEBUG} -fno-omit-frame-pointer -fsanitize=address")

add_executable(DAF main.cc ${SOURCES})