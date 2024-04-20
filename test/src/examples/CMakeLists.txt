file(GLOB_RECURSE XFDTD_EXAMPLE_SOURCES "${CMAKE_CURRENT_SOURCE_DIR}/*/*.cpp" "${CMAKE_CURRENT_SOURCE_DIR}/*.cpp")

add_custom_target(build-examples COMMAND ${CMAKE_COMMAND} -E echo "Build all examples")

set(XFDTD_EXAMPLES_LIST "")

# make example xxx
foreach(xfdtd_example_source ${XFDTD_EXAMPLE_SOURCES})
    get_filename_component(xfdtd_example_file_name ${xfdtd_example_source} NAME)
    string(REPLACE ".cpp" "" xfdtd_example_name ${xfdtd_example_file_name})
    add_executable(${xfdtd_example_name} EXCLUDE_FROM_ALL ${xfdtd_example_source})
    target_link_libraries(${xfdtd_example_name} xfdtd_core)

    add_test(NAME ${xfdtd_example_name}_test COMMAND ${xfdtd_example_name})

    add_dependencies(build-examples ${xfdtd_example_name})
endforeach()
