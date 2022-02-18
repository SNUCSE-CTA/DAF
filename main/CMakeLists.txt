add_compile_options(-O3)
add_compile_options(-Wall)
add_compile_options(-std=c++17)
add_compile_options(-DNDEBUG)

# add_compile_options(-fsanitize=address)
# add_link_options(-fsanitize=address)

add_executable(DAF main.cc ${SOURCES})
