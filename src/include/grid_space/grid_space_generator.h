#ifndef _XFDTD_LIB_GRID_SPACE_GENERATOR_H_
#define _XFDTD_LIB_GRID_SPACE_GENERATOR_H_

#include <xfdtd/boundary/boundary.h>
#include <xfdtd/grid_space/grid_space.h>

#include <memory>

namespace xfdtd {

class GridSpaceGenerator {
 public:
  static std::unique_ptr<GridSpace> generate(
      const std::vector<std::shared_ptr<Shape>>& shapes, double based_dx,
      double based_dy, double based_dz);

  static std::unique_ptr<GridSpace> generate(
      const std::vector<std::shared_ptr<Shape>>& shapes,
      const std::vector<std::shared_ptr<Boundary>>& boundaries, double based_dx,
      double based_dy, double based_dz);

 private:
  static std::unique_ptr<GridSpace> generateGridSpace1D(
      const std::vector<std::shared_ptr<Shape>>& shapes, double dz);

  static std::unique_ptr<GridSpace> generateGridSpace2D(
      const std::vector<std::shared_ptr<Shape>>& shapes, double dx, double dy);

  static std::unique_ptr<GridSpace> generateGridSpace3D(
      const std::vector<std::shared_ptr<Shape>>& shapes, double dx, double dy,
      double dz);

  static GridSpace::Dimension decideDimension(
      const std::vector<std::shared_ptr<Shape>>& shapes);
};

}  // namespace xfdtd

#endif  // _XFDTD_LIB_GRID_SPACE_GENERATOR_H_
