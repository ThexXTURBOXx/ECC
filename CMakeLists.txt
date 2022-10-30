cmake_minimum_required(VERSION 3.18)
project(ECC LANGUAGES CXX)

set(CMAKE_CXX_STANDARD 14)

find_library(COCOA_LIBRARY NAMES CoCoA cocoa REQUIRED)

add_library(util src/util/utils.C)
add_executable(BCH src/BCH/test.C)

include_directories(include)

target_link_libraries(BCH PRIVATE ${COCOA_LIBRARY} util -lgmp)
