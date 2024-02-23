#include <xfdtd/object/lumped_element/pec_plane.h>

#include <xtensor.hpp>

#include "xfdtd/object/object.h"

namespace xfdtd {

PecPlane::PecPlane(std::string name, std::unique_ptr<Cube> cube)
    : Object{std::move(name), std::move(cube), Material::createPec()} {}

}  // namespace xfdtd
