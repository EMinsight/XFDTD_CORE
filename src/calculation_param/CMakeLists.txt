aux_source_directory(. XFDTD_CALCULATION_PARAM_SRC)
add_library(xfdtd_calculation_param OBJECT ${XFDTD_CALCULATION_PARAM_SRC})
set(${ALL_OBJECT_FILES} ${ALL_OBJECT_FILES} $<TARGET_OBJECTS:xfdtd_calculation_param> PARENT_SCOPE)
