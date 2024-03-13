#ifndef _XFDTD_CORE_GRID_SPACE_GENERATOR_H_
#define _XFDTD_CORE_GRID_SPACE_GENERATOR_H_

#include <xfdtd/grid_space/grid_space.h>

#include <memory>
#include <vector>

namespace xfdtd {

class Boundary;

class GridSpaceGenerator {
 public:
  static std::unique_ptr<GridSpace> generateUniformGridSpace(
      const std::vector<const Shape*>& shapes, double based_dx, double based_dy,
      double based_dz);

 private:
  static std::unique_ptr<GridSpace> generateGridSpace1D(const Cube* domain,
                                                        double dz);

  static std::unique_ptr<GridSpace> generateGridSpace2D(const Cube* domain,
                                                        double dx, double dy);

  static std::unique_ptr<GridSpace> generateGridSpace3D(const Cube* domain,
                                                        double dx, double dy,
                                                        double dz);

  static GridSpace::Dimension decideDimension(const Shape* shape);
};

}  // namespace xfdtd

#endif  // _XFDTD_CORE_GRID_SPACE_GENERATOR_H_
