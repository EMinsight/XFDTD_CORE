aux_source_directory(. XFDTD_NFFFT_SRC)
add_library(xfdtd_nffft OBJECT ${XFDTD_NFFFT_SRC})
set(${ALL_OBJECT_FILES} ${ALL_OBJECT_FILES} $<TARGET_OBJECTS:xfdtd_nffft> PARENT_SCOPE)
