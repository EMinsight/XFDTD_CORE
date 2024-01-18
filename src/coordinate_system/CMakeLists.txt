aux_source_directory(. XFDTD_COORDINATE_SYSTEM_SRC)
add_library(xfdtd_coordinate_system OBJECT ${XFDTD_COORDINATE_SYSTEM_SRC})
set(${ALL_OBJECT_FILES} ${ALL_OBJECT_FILES} $<TARGET_OBJECTS:xfdtd_coordinate_system> PARENT_SCOPE)
