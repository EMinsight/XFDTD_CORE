aux_source_directory(. XFDTD_WAVEFORM_SRC)
add_library(xfdtd_waveform OBJECT ${XFDTD_WAVEFORM_SRC})
set(${ALL_OBJECT_FILES} ${ALL_OBJECT_FILES} $<TARGET_OBJECTS:xfdtd_waveform> PARENT_SCOPE)
