aux_source_directory(. XFDTD_SIMULATION_SRC)
add_library(xfdtd_simulation OBJECT ${XFDTD_SIMULATION_SRC})
set(${ALL_OBJECT_FILES} ${ALL_OBJECT_FILES} $<TARGET_OBJECTS:xfdtd_simulation> PARENT_SCOPE)
