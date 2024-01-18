aux_source_directory(. XFDTD_LUMPED_ELEMENT_SRC)
add_library(xfdtd_lumped_element OBJECT ${XFDTD_LUMPED_ELEMENT_SRC})
set(${ALL_OBJECT_FILES} ${ALL_OBJECT_FILES} $<TARGET_OBJECTS:xfdtd_lumped_element> PARENT_SCOPE)
