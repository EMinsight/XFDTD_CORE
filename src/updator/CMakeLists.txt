aux_source_directory(. XFDTD_UPDATOR_SRC)
aux_source_directory(${CMAKE_CURRENT_SOURCE_DIR}/ade_updator XFDTD_ADE_UPDATOR_SCR)
add_library(xfdtd_updator OBJECT ${XFDTD_UPDATOR_SRC} ${XFDTD_ADE_UPDATOR_SCR})
set(${ALL_OBJECT_FILES} ${ALL_OBJECT_FILES} $<TARGET_OBJECTS:updator> PARENT_SCOPE)
