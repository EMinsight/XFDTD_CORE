#include <xfdtd/shape/cube.h>
#include <xfdtd/shape/shape.h>

namespace xfdtd {

std::string Shape::toString() const { return std::string{"Shape()"}; }

std::unique_ptr<Cube> Shape::makeWrappedCube(
    const std::vector<const Shape*>& shapes) {
  if (shapes.empty()) {
    throw XFDTDShapeException{"Shape::makeWrappedCube: Shapes are empty"};
  }

  auto min_x = std::numeric_limits<double>::max();
  auto max_x = std::numeric_limits<double>::min();
  auto min_y = std::numeric_limits<double>::max();
  auto max_y = std::numeric_limits<double>::min();
  auto min_z = std::numeric_limits<double>::max();
  auto max_z = std::numeric_limits<double>::min();

  for (const auto& s : shapes) {
    if (s == nullptr) {
      throw XFDTDShapeException{"Shape::makeWrappedCube: Shape is nullptr"};
    }

    auto cube = s->wrappedCube();
    min_x = std::min(min_x, cube->originX());
    max_x = std::max(max_x, cube->endX());
    min_y = std::min(min_y, cube->originY());
    max_y = std::max(max_y, cube->endY());
    min_z = std::min(min_z, cube->originZ());
    max_z = std::max(max_z, cube->endZ());
  }

  return std::make_unique<Cube>(
      Vector{min_x, min_y, min_z},
      Vector{max_x - min_x, max_y - min_y, max_z - min_z});
}

}  // namespace xfdtd
