cmake_minimum_required(VERSION 3.20)

project(XFDTD_LIB_TEST VERSION 0.0.0 LANGUAGES CXX)

add_subdirectory(src)

add_test(NAME COORDINATE_SYSTEM_TEST COMMAND coordinate_system_test)
add_test(NAME SHAPE_TEST COMMAND shape_test)
