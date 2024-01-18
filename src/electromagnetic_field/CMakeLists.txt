aux_source_directory(. XFDTD_EMF_SRC)
add_library(xfdtd_emf OBJECT ${XFDTD_EMF_SRC})
set(${ALL_OBJECT_FILES} ${ALL_OBJECT_FILES} $<TARGET_OBJECTS:xfdtd_emf> PARENT_SCOPE)
