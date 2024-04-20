cmake_minimum_required(VERSION 3.20)

project(
  XFDTD_CORE
  VERSION 0.0.0
  LANGUAGES CXX)

set(XFDTD_CORE_MAIN_PROJECT OFF)
if(CMAKE_SOURCE_DIR STREQUAL PROJECT_SOURCE_DIR)
  set(XFDTDLib_MAIN_PROJECT ON)
endif()

if(PROJECT_SOURCE_DIR STREQUAL PROJECT_BINARY_DIR)
  message(FATAL_ERROR "In-source builds are not allowed")
endif()

set(CMAKE_CXX_STANDARD 20)
set(CMAKE_CXX_STANDARD_REQUIRED ON)
set(CMAKE_CXX_EXTENSIONS OFF)
set(CMAKE_EXPORT_COMPILE_COMMANDS ON)

if(NOT CMAKE_BUILD_TYPE)
  message(STATUS "No build type selected, default to Release")
  set(CMAKE_BUILD_TYPE "Release")
endif()

# option to enable MPI support commnad: cmake -D XFDTD_CORE_WITH_MPI=ON
option(XFDTD_CORE_WITH_MPI "Enable MPI support" OFF)
# option to enable single precision commnad: cmake -D
# XFDTD_CORE_SINGLE_PRECISION=ON
option(XFDTD_CORE_SINGLE_PRECISION "Enable single precision" OFF)
if(XFDTD_CORE_SINGLE_PRECISION)
  message(STATUS "XFDTD Core Single precision is enabled")
  add_definitions(-DXFDTD_CORE_SINGLE_PRECISION)
endif()

# ##############################################################################
# DEPENDENCIES
# ##############################################################################

if(WIN32)
  find_package(OpenBLAS REQUIRED)
  set(BLAS_LIBRARIES ${CMAKE_INSTALL_PREFIX} ${OpenBLAS_LIBRARIES})
else()
  find_package(BLAS REQUIRED)
  find_package(LAPACK REQUIRED)
endif()

find_package(xtl REQUIRED)
find_package(xtensor REQUIRED)
find_package(xtensor-blas REQUIRED)

set(XFDTD_CORE_PUBLIC_INCLUDE_DIR ${PROJECT_SOURCE_DIR}/include)
set(XFDTD_CORE_PRIVATE_INCLUDE_DIR ${PROJECT_SOURCE_DIR}/src/include)
set(XFDTD_CORE_INSTALL_PREFIX ${CMAKE_INSTALL_PREFIX})

include_directories(
  ${XFDTD_CORE_PUBLIC_INCLUDE_DIR} ${XFDTD_CORE_PRIVATE_INCLUDE_DIR}
  ${xtensor_INCLUDE_DIRS} ${xtensor_blas_INCLUDE_DIRS})

set(XFDTD_CORE_THIRD_PARTY_LIB ${BLAS_LIBRARIES} ${LAPACK_LIBRARIES} xtensor)

if(XFDTD_CORE_WITH_MPI)
  message(STATUS "XFDTD Core MPI support enabled")
  find_package(MPI REQUIRED)
  add_definitions(-DXFDTD_CORE_WITH_MPI)
  set(XFDTD_CORE_THIRD_PARTY_LIB ${XFDTD_CORE_THIRD_PARTY_LIB} MPI::MPI_CXX)
endif()

enable_testing()

# Output directory.
set(CMAKE_ARCHIVE_OUTPUT_DIRECTORY ${CMAKE_BINARY_DIR}/lib)
set(CMAKE_LIBRARY_OUTPUT_DIRECTORY ${CMAKE_BINARY_DIR}/lib)
set(CMAKE_RUNTIME_OUTPUT_DIRECTORY ${CMAKE_BINARY_DIR}/bin)

add_subdirectory(src)
add_subdirectory(test)

install(
  EXPORT xfdtd_core-targets
  FILE xfdtd_core-targets.cmake
  NAMESPACE xfdtd::
  DESTINATION ${XFDTD_CORE_INSTALL_PREFIX}/lib/cmake/xfdtd_core)

include(CMakePackageConfigHelpers)

configure_package_config_file(
  cmake/xfdtd_core-config.cmake.in
  ${CMAKE_CURRENT_BINARY_DIR}/xfdtd_core-config.cmake
  INSTALL_DESTINATION ${XFDTD_CORE_INSTALL_PREFIX}/lib/cmake/xfdtd_core)

write_basic_package_version_file(
  ${CMAKE_CURRENT_BINARY_DIR}/xfdtd_core-config_version.cmake
  VERSION ${PROJECT_VERSION}
  COMPATIBILITY SameMajorVersion)

install(FILES ${CMAKE_CURRENT_BINARY_DIR}/xfdtd_core-config.cmake
              ${CMAKE_CURRENT_BINARY_DIR}/xfdtd_core-config_version.cmake
        DESTINATION ${XFDTD_CORE_INSTALL_PREFIX}/lib/cmake/xfdtd_core)
