cmake_minimum_required(VERSION 2.8.9)

# Maybe stop from CMAKEing in the wrong place
if (CMAKE_BINARY_DIR STREQUAL CMAKE_SOURCE_DIR)
    message(FATAL_ERROR "Source and build directories cannot be the same. Go use the /build directory.")
endif()

include_directories(include)

# Create a single library for the project
add_library(
        rectangle-bin-pack
        src/GuillotineBinPack.cpp
        src/MaxRectsBinPack.cpp
        src/Rect.cpp
        src/ShelfBinPack.cpp
        src/ShelfNextFitBinPack.cpp
        src/SkylineBinPack.cpp
        )