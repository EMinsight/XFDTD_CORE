cmake_minimum_required(VERSION 3.20)

project(
  XFDTD_CORE
  VERSION 0.1.0
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

option(XFDTD_CORE_WITH_MPI "Enable MPI support" OFF)
# XFDTD_CORE_SINGLE_PRECISION=ON
option(XFDTD_CORE_SINGLE_PRECISION "Enable single precision" OFF)

# ##############################################################################
# DEPENDENCIES
# ##############################################################################

find_package(xtensor REQUIRED)
set(XFDTD_CORE_REQUIRED_MAX_XTENSOR_VERSION 0.25.0)
if (xtensor_FOUND)
  if (xtensor_VERSION VERSION_GREATER XFDTD_CORE_REQUIRED_MAX_XTENSOR_VERSION)
    message(FATAL_ERROR "xtensor version ${xtensor_VERSION} is not supported, please use version < ${XFDTD_CORE_REQUIRED_MAX_XTENSOR_VERSION}")
  endif()
endif()

set(XFDTD_CORE_PUBLIC_DEFINATIONS "")
set(XFDTD_CORE_PUBLIC_INCLUDE_DIR ${PROJECT_SOURCE_DIR}/include)
set(XFDTD_CORE_PRIVATE_INCLUDE_DIR ${PROJECT_SOURCE_DIR}/src/include)
set(XFDTD_CORE_THIRD_PARTY_INCLUDE_DIR ${xtensor_INCLUDE_DIRS})
set(XFDTD_CORE_INSTALL_PREFIX ${CMAKE_INSTALL_PREFIX})

set(XFDTD_CORE_THIRD_PARTY_LIB ${XFDTD_CORE_THIRD_PARTY_LIB} xtensor)

set(XFDTD_CORE_NEED_LINKER_TBB OFF)
if(CMAKE_CXX_COMPILER_ID MATCHES "MSVC")
  # VS 2017 or later
  set(XFDTD_CORE_DEFINATIONS ${XFDTD_CORE_DEFINATIONS} XFDTD_CORE_PSTL_ENABLE)
  message(STATUS "XFDTD Core PSTL support is enabled")
elseif(CMAKE_CXX_COMPILER_ID MATCHES "GNU" OR CMAKE_CXX_COMPILER_ID MATCHES "Clang")
  if(CMAKE_SYSTEM_NAME MATCHES "Windows" OR CMAKE_SYSTEM_NAME MATCHES "Darwin")
    # TODO: mingw64 works
    set(XFDTD_CORE_NEED_LINKER_TBB OFF)
  else()
    set(XFDTD_CORE_NEED_LINKER_TBB ON)
  endif()
endif()

if (XFDTD_CORE_NEED_LINKER_TBB)
  find_package(TBB QUIET)
  if (NOT TBB_FOUND)
    set(XFDTD_CORE_NEED_LINKER_TBB OFF)
  else()
    set(XFDTD_CORE_DEFINATIONS ${XFDTD_CORE_DEFINATIONS} XFDTD_CORE_PSTL_ENABLE)
    set(XFDTD_CORE_THIRD_PARTY_LIB ${XFDTD_CORE_THIRD_PARTY_LIB} TBB::tbb)
    message(STATUS "XFDTD Core PSTL support is enabled")
  endif()
endif()

if(XFDTD_CORE_WITH_MPI)
  message(STATUS "XFDTD Core MPI support is enabled")
  find_package(MPI REQUIRED)
  set(XFDTD_CORE_DEFINATIONS ${XFDTD_CORE_DEFINATIONS} XFDTD_CORE_WITH_MPI)
  set(XFDTD_CORE_THIRD_PARTY_LIB ${XFDTD_CORE_THIRD_PARTY_LIB} MPI::MPI_CXX)
  set(XFDTD_CORE_THIRD_PARTY_INCLUDE_DIR ${XFDTD_CORE_THIRD_PARTY_INCLUDE_DIR}
                                         ${MPI_CXX_INCLUDE_DIRS})
endif()
if(XFDTD_CORE_SINGLE_PRECISION)
  message(STATUS "XFDTD Core Single precision is enabled")
  set(XFDTD_CORE_DEFINATIONS ${XFDTD_CORE_DEFINATIONS} XFDTD_CORE_SINGLE_PRECISION)
endif()

add_compile_definitions(${XFDTD_CORE_DEFINATIONS})
include_directories(
  ${XFDTD_CORE_PUBLIC_INCLUDE_DIR} ${XFDTD_CORE_PRIVATE_INCLUDE_DIR}
  ${XFDTD_CORE_THIRD_PARTY_INCLUDE_DIR})

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
