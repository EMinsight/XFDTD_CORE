#include "xfdtd/material/dispersive_material.h"

namespace xfdtd {

LinearDispersiveMaterial::LinearDispersiveMaterial(std::string_view name,
                                                   Type type,
                                                   ElectroMagneticProperty emp)
    : Material{name, emp, true}, _type{type} {}

}  // namespace xfdtd
