cmake_minimum_required(VERSION 2.8.9)

# Create a library
add_library(tinyply src/tinyply.cpp)
target_include_directories(tinyply PUBLIC "${CMAKE_CURRENT_SOURCE_DIR}/include/")
