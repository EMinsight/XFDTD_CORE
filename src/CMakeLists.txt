add_subdirectory(coordinate_system)
add_subdirectory(shape)
add_subdirectory(grid_space)
add_subdirectory(calculation_param)
add_subdirectory(electromagnetic_field)
add_subdirectory(material)
add_subdirectory(object)
add_subdirectory(updator)
add_subdirectory(waveform)
add_subdirectory(waveform_source)
add_subdirectory(boundary)
add_subdirectory(monitor)
add_subdirectory(network)
add_subdirectory(nffft)
add_subdirectory(simulation)
add_subdirectory(domain)

add_library(xfdtd STATIC ${ALL_OBJECT_FILES})

set(XFDTD_LIBS
    xfdtd_coordinate_system
    xfdtd_shape
    xfdtd_grid_space
    xfdtd_calculation_param
    xfdtd_emf
    xfdtd_material
    xfdtd_object
    xfdtd_lumped_element
    xfdtd_updator
    xfdtd_waveform
    xfdtd_waveform_source
    xfdtd_boundary
    xfdtd_monitor
    xfdtd_network
    xfdtd_nffft
    xfdtd_simulation
    xfdtd_domain)

target_link_libraries(
  xfdtd PUBLIC ${XFDTD_LIBS} ${XFDTD_LIB_THIRD_PARTY_LIB})

target_include_directories(
  xfdtd PUBLIC $<BUILD_INTERFACE:${XFDTD_LIB_PUBLIC_INCLUDE_DIR}>
               $<INSTALL_INTERFACE:include/xfdtd>)

target_include_directories(
  xfdtd PRIVATE $<BUILD_INTERFACE:${XFDTD_LIB_PRIVATE_INCLUDE_DIR}>
                $<INSTALL_INTERFACE:include>)

install(
  TARGETS xfdtd ${XFDTD_LIBS}
  EXPORT xfdtd-targets
  ARCHIVE DESTINATION ${XFDTD_INSTALL_PREFIX}/lib
  LIBRARY DESTINATION ${XFDTD_INSTALL_PREFIX}/lib
  RUNTIME DESTINATION ${XFDTD_INSTALL_PREFIX}/bin
  INCLUDES
  DESTINATION ${XFDTD_INSTALL_PREFIX}/include)

install(DIRECTORY ${XFDTD_LIB_PUBLIC_INCLUDE_DIR}/xfdtd DESTINATION ${XFDTD_INSTALL_PREFIX}/include)
