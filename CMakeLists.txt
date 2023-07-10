cmake_minimum_required(VERSION 3.18)
project(ECC LANGUAGES CXX)

set(CMAKE_CXX_STANDARD 14)

find_library(COCOA_LIBRARY NAMES CoCoA cocoa REQUIRED)

# Util libraries
add_library(util src/util/utils.C)
target_link_libraries(util PRIVATE ${COCOA_LIBRARY} -lgmp)

# ECC types
add_library(types src/types/cyclic.C src/types/linear.C)
target_link_libraries(types PRIVATE ${COCOA_LIBRARY} util -lgmp)

# ECCs
add_library(BCH src/BCH/bch.C)
target_link_libraries(BCH PRIVATE ${COCOA_LIBRARY} util types -lgmp)
add_library(Golay src/Golay/golay.C)
target_link_libraries(Golay PRIVATE ${COCOA_LIBRARY} util types -lgmp)
add_library(Ham src/Ham/ham.C)
target_link_libraries(Ham PRIVATE ${COCOA_LIBRARY} util types -lgmp)
add_library(RM src/RM/rm.C)
target_link_libraries(RM PRIVATE ${COCOA_LIBRARY} util types -lgmp)

# Fuzzy Extractors
add_library(Fuzzy src/fuzzy/fuzzy.C SHA256/src/SHA256.cpp)
target_link_libraries(Fuzzy PRIVATE ${COCOA_LIBRARY} util types -lgmp)

# Executables
include_directories(include SHA256/include)
add_executable(test src/test/test.C src/test/ecc_test.C src/test/fuzzy_test.C)
target_link_libraries(test PRIVATE ${COCOA_LIBRARY} util types BCH Golay Ham RM Fuzzy -lgmp)
