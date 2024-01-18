aux_source_directory(. XFDTD_SHAPE_SRC)
add_library(xfdtd_shape OBJECT ${XFDTD_SHAPE_SRC})
set(${ALL_OBJECT_FILES} ${ALL_OBJECT_FILES} $<TARGET_OBJECTS:xfdtd_shape> PARENT_SCOPE)
