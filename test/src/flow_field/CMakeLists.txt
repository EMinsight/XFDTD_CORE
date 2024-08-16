aux_source_directory(. FLOW_FIELD_SRC)
message(STATUS "FLOW_FIELD_SRC: ${FLOW_FIELD_SRC}")
add_executable(flow_field ${FLOW_FIELD_SRC})
target_link_libraries(flow_field PRIVATE xfdtd_core)