aux_source_directory(. XFDTD_WAVEFORM_SOURCE_SRC)
add_library(xfdtd_waveform_source OBJECT ${XFDTD_WAVEFORM_SOURCE_SRC})
set(${ALL_OBJECT_FILES} ${ALL_OBJECT_FILES} $<TARGET_OBJECTS:xfdtd_waveform_source> PARENT_SCOPE)
