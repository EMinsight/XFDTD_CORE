#include <xfdtd/object/flow_field.h>

#include <memory>

#include "xfdtd/common/type_define.h"
#include "xfdtd/coordinate_system/coordinate_system.h"
#include "xfdtd/shape/cube.h"

namespace xfdtd {

auto FlowField::correctUpdateCoefficient() -> void {
  const auto& grid_space = gridSpace();

  if (grid_space->type() != GridSpace::Type::UNIFORM) {
    throw XFDTDFlowFieldException{
        "handleDispersion(): Non-uniform grid space is not supported yet"};
  }

  const auto& flow_field_entries = _flow_field_entries;
  const auto& region =
      Cube{Vector{grid_space->eNodeX().front(), grid_space->eNodeY().front(),
                  grid_space->eNodeZ().front()},
           Vector{grid_space->eNodeX().back() - grid_space->eNodeX().front(),
                  grid_space->eNodeY().back() - grid_space->eNodeY().front(),
                  grid_space->eNodeZ().back() - grid_space->eNodeZ().front()}};

  for (const auto& e : flow_field_entries) {
    auto vec = e.position();
    if (!region.isInside(vec, 1e-6)) {
      continue;
    }

    auto g = grid_space->getGrid(vec);
    auto i = g.i();
    auto j = g.j();
    auto k = g.k();

    correctCoefficient(i, j, k, e);
  }
}

auto FlowField::correctCoefficient(Index i, Index j, Index k,
                                   const FlowFieldEntry& entry) -> void {}

}  // namespace xfdtd