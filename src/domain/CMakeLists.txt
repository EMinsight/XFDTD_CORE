aux_source_directory(. XFDTD_DOMAIN_SRC)
add_library(xfdtd_domain OBJECT ${XFDTD_DOMAIN_SRC})
set(${ALL_OBJECT_FILES} ${ALL_OBJECT_FILES} $<TARGET_OBJECTS:xfdtd_domain PARENT_SCOPE)
