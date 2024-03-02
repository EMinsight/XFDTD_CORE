cmake_minimum_required(VERSION 3.20)

project(XFDTDLib VERSION 0.0.0 LANGUAGES CXX)

set(XFDTDLib_MAIN_PROJECT OFF)
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

set(XFDTD_LIB_PUBLIC_INCLUDE_DIR ${PROJECT_SOURCE_DIR}/include)
set(XFDTD_LIB_PRIVATE_INCLUDE_DIR ${PROJECT_SOURCE_DIR}/src/include)
set(XFDTD_INSTALL_PREFIX ${CMAKE_INSTALL_PREFIX})

include_directories(${XFDTD_LIB_PUBLIC_INCLUDE_DIR}
                    ${XFDTD_LIB_PRIVATE_INCLUDE_DIR}
                    ${FFTW_INCLUDE_DIR}
                    ${xtensor_INCLUDE_DIRS}
                    ${xtensor_blas_INCLUDE_DIRS}
                    ${xtensor_fftw_INCLUDE_DIRS})

set(XFDTD_LIB_THIRD_PARTY_LIB ${FFTW_LIBRARIES}
                              ${BLAS_LIBRARIES}
                              ${LAPACK_LIBRARIES}
                              xtensor
                              xtensor-blas)

enable_testing()

# # for all target add -fsanitize=thread
# if (XFDTDLib_MAIN_PROJECT)
#   set(CMAKE_CXX_FLAGS "${CMAKE_CXX_FLAGS} -fsanitize=thread")
#   set(CMAKE_CXX_FLAGS_DEBUG "${CMAKE_CXX_FLAGS_DEBUG} -fsanitize=thread")
#   set(CMAKE_LINKER_FLAGS "${CMAKE_LINKER_FLAGS} -fsanitize=thread")
#   set(CMAKE_LINKER_FLAGS_DEBUG "${CMAKE_LINKER_FLAGS_DEBUG} -fsanitize=thread")
# endif()

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
  cmake/xfdtd-config.cmake.in
  ${CMAKE_CURRENT_BINARY_DIR}/xfdtd-config.cmake
  INSTALL_DESTINATION ${XFDTD_INSTALL_PREFIX}/lib/cmake/xfdtd)

write_basic_package_version_file(
  ${CMAKE_CURRENT_BINARY_DIR}/xfdtd-config-version.cmake
  VERSION ${PROJECT_VERSION}
  COMPATIBILITY SameMajorVersion)

install(FILES
  ${CMAKE_CURRENT_BINARY_DIR}/xfdtd-config.cmake
  ${CMAKE_CURRENT_BINARY_DIR}/xfdtd-config-version.cmake
  DESTINATION ${XFDTD_INSTALL_PREFIX}/lib/cmake/xfdtd
)
