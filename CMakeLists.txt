cmake_minimum_required(VERSION 3.18)
project(ECC LANGUAGES CXX)

find_library(COCOA_LIBRARY NAMES CoCoA cocoa REQUIRED)

add_executable(RS src/RS/test.C)

target_link_libraries(RS PRIVATE ${COCOA_LIBRARY} -lgmp)
