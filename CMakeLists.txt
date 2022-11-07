cmake_minimum_required(VERSION 3.18)
project(ECC LANGUAGES CXX)

set(CMAKE_CXX_STANDARD 14)

find_library(COCOA_LIBRARY NAMES CoCoA cocoa REQUIRED)

# Util libraries
add_library(util src/util/utils.C)

# ECC types
add_library(types src/types/cyclic.C src/types/linear.C)

# ECCs
add_library(BCH src/BCH/bch.C)
target_link_libraries(BCH PRIVATE ${COCOA_LIBRARY} util types -lgmp)
add_library(Ham src/Ham/ham.C)
target_link_libraries(Ham PRIVATE ${COCOA_LIBRARY} util types -lgmp)
add_library(RM src/RM/rm.C)
target_link_libraries(RM PRIVATE ${COCOA_LIBRARY} util types -lgmp)

# Executables
add_executable(test src/test/test.C)

include_directories(include)

target_link_libraries(test PRIVATE ${COCOA_LIBRARY} util types BCH Ham RM -lgmp)
