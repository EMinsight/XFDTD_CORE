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

# option to enable MPI support 
# commnad: cmake -D XFDTD_CORE_WITH_MPI=ON ${PROJECT_BIUILD_DIR}
option(XFDTD_CORE_WITH_MPI "Enable MPI support" OFF)

include("cmake/FindFFTW.cmake")
if(WIN32)
  find_package(OpenBLAS REQUIRED)
  set(BLAS_LIBRARIES ${CMAKE_INSTALL_PREFIX}${OpenBLAS_LIBRARIES})
else()
  find_package(BLAS REQUIRED)
  find_package(LAPACK REQUIRED)
endif()

find_package(xtl REQUIRED)
find_package(xtensor REQUIRED)
find_package(xtensor-blas REQUIRED)
set(FFTW_USE_DOUBLE ON)
find_package(xtensor-fftw REQUIRED)

set(XFDTD_CORE_PUBLIC_INCLUDE_DIR ${PROJECT_SOURCE_DIR}/include)
set(XFDTD_CORE_PRIVATE_INCLUDE_DIR ${PROJECT_SOURCE_DIR}/src/include)
set(XFDTD_CORE_INSTALL_PREFIX ${CMAKE_INSTALL_PREFIX})

include_directories(
  ${XFDTD_CORE_PUBLIC_INCLUDE_DIR} ${XFDTD_CORE_PRIVATE_INCLUDE_DIR}
  ${FFTW_INCLUDE_DIR} ${xtensor_INCLUDE_DIRS} ${xtensor_blas_INCLUDE_DIRS}
  ${xtensor_fftw_INCLUDE_DIRS})

set(XFDTD_CORE_THIRD_PARTY_LIB ${FFTW_LIBRARIES} ${BLAS_LIBRARIES}
                               ${LAPACK_LIBRARIES} xtensor)

if(XFDTD_CORE_WITH_MPI)
  message(STATUS "XFDTD MPI support enabled")
  find_package(MPI REQUIRED)
  add_definitions(-DXFDTD_CORE_WITH_MPI)
  include_directories(${MPI_CXX_INCLUDE_PATH})
  set(XFDTD_CORE_THIRD_PARTY_LIB ${XFDTD_CORE_THIRD_PARTY_LIB}
                                 ${MPI_CXX_LIBRARIES})
endif()

enable_testing()

add_subdirectory(src)
add_subdirectory(test)
add_subdirectory(examples)

install(
  EXPORT xfdtd-targets
  FILE xfdtd-targets.cmake
  NAMESPACE xfdtd::
  DESTINATION ${XFDTD_INSTALL_PREFIX}/lib/cmake/xfdtd)

include(CMakePackageConfigHelpers)

configure_package_config_file(
  cmake/xfdtd-config.cmake.in ${CMAKE_CURRENT_BINARY_DIR}/xfdtd-config.cmake
  INSTALL_DESTINATION ${XFDTD_INSTALL_PREFIX}/lib/cmake/xfdtd)

write_basic_package_version_file(
  ${CMAKE_CURRENT_BINARY_DIR}/xfdtd-config-version.cmake
  VERSION ${PROJECT_VERSION}
  COMPATIBILITY SameMajorVersion)

install(FILES ${CMAKE_CURRENT_BINARY_DIR}/xfdtd-config.cmake
              ${CMAKE_CURRENT_BINARY_DIR}/xfdtd-config-version.cmake
        DESTINATION ${XFDTD_INSTALL_PREFIX}/lib/cmake/xfdtd)
