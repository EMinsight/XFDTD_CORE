#include "xfdtd/material/dispersive_material.h"

namespace xfdtd {

LinearDispersiveMaterial::LinearDispersiveMaterial(const std::string& name,
                                                   Type type,
                                                   ElectroMagneticProperty emp)
    : Material{name, emp, true}, _type{type} {}

auto LinearDispersiveMaterial::type() const -> Type { return _type; }

}  // namespace xfdtd
