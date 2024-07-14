add_subdirectory(boundary)
add_subdirectory(calculation_param)
add_subdirectory(coordinate_system)
add_subdirectory(domain)
add_subdirectory(electromagnetic_field)
add_subdirectory(grid_space)
add_subdirectory(material)
add_subdirectory(monitor)
add_subdirectory(network)
add_subdirectory(nffft)
add_subdirectory(object)
add_subdirectory(parallel)
add_subdirectory(shape)
add_subdirectory(simulation)
add_subdirectory(updator)
add_subdirectory(waveform)
add_subdirectory(waveform_source)

add_library(xfdtd_core ${ALL_OBJECT_FILES})

set(XFDTD_CORE_LIBS
    xfdtd_boundary
    xfdtd_calculation_param
    xfdtd_coordinate_system
    xfdtd_domain
    xfdtd_emf
    xfdtd_grid_space
    xfdtd_lumped_element
    xfdtd_material
    xfdtd_monitor
    xfdtd_network
    xfdtd_nffft
    xfdtd_shape
    xfdtd_object
    xfdtd_parallel
    xfdtd_simulation
    xfdtd_updator
    xfdtd_waveform
    xfdtd_waveform_source
)

target_link_libraries(xfdtd_core PUBLIC ${XFDTD_CORE_THIRD_PARTY_LIB}
                                        ${XFDTD_CORE_LIBS})

target_include_directories(
  xfdtd_core PUBLIC $<BUILD_INTERFACE:${XFDTD_CORE_PUBLIC_INCLUDE_DIR}>
                    $<INSTALL_INTERFACE:include>)

target_include_directories(
  xfdtd_core PRIVATE $<BUILD_INTERFACE:${XFDTD_CORE_PRIVATE_INCLUDE_DIR}>
                     $<INSTALL_INTERFACE:include>)

install(
  TARGETS xfdtd_core ${XFDTD_CORE_LIBS}
  EXPORT xfdtd_core-targets
  ARCHIVE DESTINATION ${XFDTD_CORE_INSTALL_PREFIX}/lib
  LIBRARY DESTINATION ${XFDTD_CORE_INSTALL_PREFIX}/lib
  RUNTIME DESTINATION ${XFDTD_CORE_INSTALL_PREFIX}/bin
  INCLUDES
  DESTINATION ${XFDTD_CORE_INSTALL_PREFIX}/include)

install(DIRECTORY ${XFDTD_CORE_PUBLIC_INCLUDE_DIR}/xfdtd
        DESTINATION ${XFDTD_CORE_INSTALL_PREFIX}/include)
