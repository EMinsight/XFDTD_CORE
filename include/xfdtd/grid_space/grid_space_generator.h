#ifndef _XFDTD_CORE_GRID_SPACE_GENERATOR_H_
#define _XFDTD_CORE_GRID_SPACE_GENERATOR_H_

#include <xfdtd/common/type_define.h>
#include <xfdtd/grid_space/grid_space.h>

#include <memory>
#include <vector>

namespace xfdtd {

class Boundary;

class GridSpaceGenerator {
 public:
  static std::unique_ptr<GridSpace> generateUniformGridSpace(
      const std::vector<const Shape*>& shapes,
      const std::vector<const Boundary*>& boundaries, Real based_dx,
      Real based_dy, Real based_dz);

 private:
  static std::unique_ptr<Cube> extendDomain(
      std::unique_ptr<Cube> domain,
      const std::vector<const Boundary*>& boundaries, Real based_dx,
      Real based_dy, Real based_dz);

  static std::unique_ptr<GridSpace> generateGridSpace1D(const Cube* domain,
                                                        Real dz);

  static std::unique_ptr<GridSpace> generateGridSpace2D(const Cube* domain,
                                                        Real dx, Real dy);

  static std::unique_ptr<GridSpace> generateGridSpace3D(const Cube* domain,
                                                        Real dx, Real dy,
                                                        Real dz);

  static GridSpace::Dimension decideDimension(const Shape* shape);
};

}  // namespace xfdtd

#endif  // _XFDTD_CORE_GRID_SPACE_GENERATOR_H_
