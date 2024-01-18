aux_source_directory(. XFDTD_GRID_SPACE_SRC)
add_library(xfdtd_grid_space OBJECT ${XFDTD_GRID_SPACE_SRC})
set(${ALL_OBJECT_FILES} ${ALL_OBJECT_FILES} $<TARGET_OBJECTS:xfdtd_grid_space> PARENT_SCOPE)
