#ifndef _XFDTD_CORE_GRID_SPACE_H_
#define _XFDTD_CORE_GRID_SPACE_H_

#include <xfdtd/common/type_define.h>
#include <xfdtd/coordinate_system/coordinate_system.h>
#include <xfdtd/exception/exception.h>
#include <xfdtd/grid_space/grid.h>
#include <xfdtd/grid_space/grid_box.h>
#include <xfdtd/shape/cube.h>
#include <xfdtd/shape/shape.h>

#include <cstddef>
#include <memory>
#include <string>
#include <utility>
#include <xtensor/xindex_view.hpp>

namespace xfdtd {

template <typename Arr1D, typename Index, typename Real>
class GridSpaceData;

class XFDTDGridSpaceException : public XFDTDException {
 public:
  explicit XFDTDGridSpaceException(std::string message = "")
      : XFDTDException(std::move(message)) {}
};

class GridSpace {
 public:
  enum class Dimension { UNDEFINED, ONE, TWO, THREE };
  enum class Type { UNDEFINED, UNIFORM, NONUNIFORM };

 public:
  GridSpace(const GridSpace&) = default;

  GridSpace(GridSpace&&) noexcept = default;

  GridSpace& operator=(const GridSpace&) = default;

  GridSpace& operator=(GridSpace&&) noexcept = default;

  virtual ~GridSpace() = default;

  Dimension dimension() const;

  Type type() const;

  // Cube region() const;

  GridBox box() const;

  Real basedDx() const;

  Real basedDy() const;

  Real basedDz() const;

  Real minDx() const;

  Real minDy() const;

  Real minDz() const;

  auto& gridWithMaterial() const { return _grid_with_material; }

  std::shared_ptr<GridSpace> globalGridSpace() const;

  const Array1D<Real>& eNodeX() const;

  const Array1D<Real>& eNodeY() const;

  const Array1D<Real>& eNodeZ() const;

  const Array1D<Real>& hNodeX() const;

  const Array1D<Real>& hNodeY() const;

  const Array1D<Real>& hNodeZ() const;

  const Array1D<Real>& eSizeX() const;

  const Array1D<Real>& eSizeY() const;

  const Array1D<Real>& eSizeZ() const;

  const Array1D<Real>& hSizeX() const;

  const Array1D<Real>& hSizeY() const;

  const Array1D<Real>& hSizeZ() const;

  const Array1D<Real>& eSize(Axis::XYZ xyz) const;

  const Array1D<Real>& hSize(Axis::XYZ xyz) const;

  Index sizeX() const;

  Index sizeY() const;

  Index sizeZ() const;

  Grid getGrid(Real x, Real y, Real z) const;

  Grid getGrid(const Vector& vector) const;

  Vector getGridOriginVector(const Grid& grid) const;

  Vector getGridCenterVector(const Grid& grid) const;

  Vector getGridEndVector(const Grid& grid) const;

  GridBox getGridBox(const Shape* shape) const;

  virtual void correctGridSpace() = 0;

  void extendGridSpace(Axis::Direction direction, Index num, Real dl);

  virtual std::unique_ptr<GridSpace> subGridSpace(Index start_i, Index start_j,
                                                  Index start_k, Index end_i,
                                                  Index end_j,
                                                  Index end_k) const = 0;

  GridBox getGridBoxWithoutCheck(const Shape* shape) const;

  Grid getGridWithoutCheck(const Vector& vector) const;

  void generateMaterialGrid(Index nx, Index ny, Index nz);

  void setGlobalGridSpace(std::weak_ptr<GridSpace> global_grid_space);

  GridBox validGridBoxEx() const;

  GridBox validGridBoxEy() const;

  GridBox validGridBoxEz() const;

  GridBox globalBox() const;

  std::string toString() const;

  Grid transformNodeToGlobal(const Grid& grid) const;

  GridBox transformNodeToGlobal(const GridBox& box) const;

  auto gridSpaceData() const -> GridSpaceData<Array1D<Real>, Index, Real>;

  auto eps() const -> Real;

  Array1D<Real>& eNodeX();

  Array1D<Real>& eNodeY();

  Array1D<Real>& eNodeZ();

  Array1D<Real>& hNodeX();

  Array1D<Real>& hNodeY();

  Array1D<Real>& hNodeZ();

  Array1D<Real>& eSizeX();

  Array1D<Real>& eSizeY();

  Array1D<Real>& eSizeZ();

  Array1D<Real>& hSizeX();

  Array1D<Real>& hSizeY();

  Array1D<Real>& hSizeZ();

  Array1D<Real>& eSize(Axis::XYZ xyz);

  Array1D<Real>& hSize(Axis::XYZ xyz);

 protected:
  GridSpace(Real dx, Real dy, Real dz, Dimension dimension,
            Array1D<Real> e_node_x, Array1D<Real> e_node_y,
            Array1D<Real> e_node_z);

  // for subGridSpace
  GridSpace(Dimension dimension, Type type, GridBox global_box, Real based_dx,
            Real based_dy, Real based_dz, Real min_dx, Real min_dy, Real min_dz,
            Array1D<Real> e_node_x, Array1D<Real> e_node_y,
            Array1D<Real> e_node_z, Array1D<Real> h_node_x,
            Array1D<Real> h_node_y, Array1D<Real> h_node_z,
            Array1D<Real> e_size_x, Array1D<Real> e_size_y,
            Array1D<Real> e_size_z, Array1D<Real> h_size_x,
            Array1D<Real> h_size_y, Array1D<Real> h_size_z);

  void setMinDx(Real min_dx);

  void setMinDy(Real min_dy);

  void setMinDz(Real min_dz);

  virtual Index handleTransformX(Real x) const;

  virtual Index handleTransformY(Real y) const;

  virtual Index handleTransformZ(Real z) const;

  virtual Index handleTransformXWithoutCheck(Real x) const;

  virtual Index handleTransformYWithoutCheck(Real y) const;

  virtual Index handleTransformZWithoutCheck(Real z) const;

  void correctGridSpaceForOne(Real dl, const Array1D<Real>& e_node,
                              Array1D<Real>& h_node, Array1D<Real>& e_size,
                              Array1D<Real>& h_size);

  Real _min_x, _min_y, _min_z;
  Real _max_x, _max_y, _max_z;

 private:
  Dimension _dimension{Dimension::UNDEFINED};
  Type _type{Type::UNDEFINED};
  Real _based_dx, _based_dy, _based_dz;
  Real _min_dx, _min_dy, _min_dz;

  Array1D<Real> _e_node_x, _e_node_y, _e_node_z;
  Array1D<Real> _h_node_x, _h_node_y, _h_node_z;
  Array1D<Real> _e_size_x, _e_size_y, _e_size_z;
  Array1D<Real> _h_size_x, _h_size_y, _h_size_z;

  Array<std::shared_ptr<Grid>> _grid_with_material;

  std::weak_ptr<GridSpace> _global_grid_space;

  GridBox _global_box;

  Array1D<Real> calculateHNode(const Array1D<Real>& e_node);

  Array1D<Real> calculateESize(const Array1D<Real>& e_node);

  Array1D<Real> calculateHSize(const Array1D<Real>& h_node,
                               const Array1D<Real>& e_node, Real dl);
};

}  // namespace xfdtd

#endif  // _XFDTD_CORE_GRID_SPACE_H_
