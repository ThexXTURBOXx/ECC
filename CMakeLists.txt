cmake_minimum_required(VERSION 3.18)
project(ECC LANGUAGES CXX)

set(CMAKE_CXX_STANDARD 14)

find_library(COCOA_LIBRARY NAMES CoCoA cocoa REQUIRED)

# Util libraries
add_library(util src/util/utils.C)

# ECCs
add_library(BCH src/BCH/bch.C)

# Executables
add_executable(test src/test/test.C)

include_directories(include)

target_link_libraries(test PRIVATE ${COCOA_LIBRARY} util BCH -lgmp)
