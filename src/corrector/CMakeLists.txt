aux_source_directory(. XFDTD_CORRECTOR_SRC)
add_library(xfdtd_corrector OBJECT ${XFDTD_CORRECTOR_SRC})
set(${ALL_OBJECT_FILES} ${ALL_OBJECT_FILES} $<TARGET_OBJECTS:xfdtd_corrector> PARENT_SCOPE)
