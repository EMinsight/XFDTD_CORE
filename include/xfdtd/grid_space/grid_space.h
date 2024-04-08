#ifndef _XFDTD_CORE_GRID_SPACE_H_
#define _XFDTD_CORE_GRID_SPACE_H_

#include <xfdtd/coordinate_system/coordinate_system.h>
#include <xfdtd/exception/exception.h>
#include <xfdtd/shape/cube.h>
#include <xfdtd/shape/shape.h>
#include <xfdtd/common/type_define.h>

#include <cstddef>
#include <memory>
#include <string>
#include <utility>
#include <xtensor/xarray.hpp>
#include <xtensor/xindex_view.hpp>

namespace xfdtd {

class XFDTDGridSpaceException : public XFDTDException {
 public:
  explicit XFDTDGridSpaceException(
      std::string message = "XFDTD GridSpace Exception")
      : XFDTDException(std::move(message)) {}
};

class Grid {
 public:
  Grid() = default;

  Grid(std::size_t i, std::size_t j, std::size_t k,
       std::size_t material_index = -1);

  Grid(const Grid&) = default;

  Grid(Grid&&) noexcept = default;

  Grid& operator=(const Grid&) = default;

  Grid& operator=(Grid&&) noexcept = default;

  ~Grid() = default;

  Grid operator+(const Grid& grid) const;

  Grid operator-(const Grid& grid) const;

  std::size_t i() const;

  std::size_t j() const;

  std::size_t k() const;

  std::size_t materialIndex() const;

  void setMaterialIndex(std::size_t index);

  std::string toString() const;

 private:
  std::size_t _i{0}, _j{0}, _k{0}, _material_index{static_cast<size_t>(-1)};
};

class GridBox {
 public:
  GridBox() = default;

  GridBox(Grid origin, Grid size);

  GridBox(const GridBox&) = default;

  GridBox(GridBox&&) noexcept = default;

  GridBox& operator=(const GridBox&) = default;

  GridBox& operator=(GridBox&&) noexcept = default;

  ~GridBox() = default;

  Grid origin() const;

  Grid size() const;

  Grid center() const;

  Grid end() const;

  std::string toString() const;

 private:
  Grid _origin, _size;
};

class GridSpace {
 public:
  using GridSpaceRegion = Cube;
  enum class Dimension { UNDEFINED, ONE, TWO, THREE };
  enum class Type { UNDEFINED, UNIFORM, NONUNIFORM };

  GridSpace(const GridSpace&) = default;

  GridSpace(GridSpace&&) noexcept = default;

  GridSpace& operator=(const GridSpace&) = default;

  GridSpace& operator=(GridSpace&&) noexcept = default;

  virtual ~GridSpace() = default;

  Dimension dimension() const;

  Type type() const;

  GridSpaceRegion region() const;

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

  std::size_t sizeX() const;

  std::size_t sizeY() const;

  std::size_t sizeZ() const;

  Grid getGrid(Real x, Real y, Real z) const;

  Grid getGrid(const Vector& vector) const;

  Vector getGridOriginVector(const Grid& grid) const;

  Vector getGridCenterVector(const Grid& grid) const;

  Vector getGridEndVector(const Grid& grid) const;

  auto getShapeMask(const Shape* shape) const;

  template<typename T>
  auto getGridWithMaterialView(const T& mask);

  GridBox getGridBox(const Shape* shape) const;

  virtual void correctGridSpace() = 0;

  void extendGridSpace(Axis::Direction direction, std::size_t num, Real dl);

  virtual std::unique_ptr<GridSpace> subGridSpace(
      std::size_t start_i, std::size_t start_j, std::size_t start_k,
      std::size_t end_i, std::size_t end_j, std::size_t end_k) const = 0;

  GridBox getGridBoxWithoutCheck(const Shape* shape) const;

  Grid getGridWithoutCheck(const Vector& vector) const;

  void generateMaterialGrid(std::size_t nx, std::size_t ny, std::size_t nz);

  void setGlobalGridSpace(std::weak_ptr<GridSpace> global_grid_space);

  GridBox validGridBoxEx() const;

  GridBox validGridBoxEy() const;

  GridBox validGridBoxEz() const;

  GridBox globalBox() const;

  std::string toString() const;

  Grid transformNodeToGlobal(const Grid& grid) const;

  GridBox transformNodeToGlobal(const GridBox& box) const;

 protected:
  GridSpace(Real dx, Real dy, Real dz, Dimension dimension,
            Array1D<Real> e_node_x, Array1D<Real> e_node_y,
            Array1D<Real> e_node_z);

  // for subGridSpace
  GridSpace(Dimension dimension, Type type, GridBox global_box, Real based_dx,
            Real based_dy, Real based_dz, Real min_dx, Real min_dy,
            Real min_dz, Array1D<Real> e_node_x,
            Array1D<Real> e_node_y, Array1D<Real> e_node_z,
            Array1D<Real> h_node_x, Array1D<Real> h_node_y,
            Array1D<Real> h_node_z, Array1D<Real> e_size_x,
            Array1D<Real> e_size_y, Array1D<Real> e_size_z,
            Array1D<Real> h_size_x, Array1D<Real> h_size_y,
            Array1D<Real> h_size_z);

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

  void setMinDx(Real min_dx);

  void setMinDy(Real min_dy);

  void setMinDz(Real min_dz);

  virtual std::size_t handleTransformX(Real x) const;

  virtual std::size_t handleTransformY(Real y) const;

  virtual std::size_t handleTransformZ(Real z) const;

  virtual std::size_t handleTransformXWithoutCheck(Real x) const;

  virtual std::size_t handleTransformYWithoutCheck(Real y) const;

  virtual std::size_t handleTransformZWithoutCheck(Real z) const;

  void correctGridSpaceForOne(Real dl, const Array1D<Real>& e_node,
                              Array1D<Real>& h_node,
                              Array1D<Real>& e_size,
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
                                    const Array1D<Real>& e_node,
                                    Real dl);
};

inline auto GridSpace::getShapeMask(const Shape* shape) const {
  xt::xarray<bool> mask{xt::make_lambda_xfunction(
      [shape, this](auto&& g) {
        return shape->isInside(getGridCenterVector(*g));
      },
      _grid_with_material)};
  return mask;
}

template<typename T>
inline auto GridSpace::getGridWithMaterialView(const T& mask) {
  return xt::filter(_grid_with_material, mask);
}

}  // namespace xfdtd

#endif  // _XFDTD_CORE_GRID_SPACE_H_
