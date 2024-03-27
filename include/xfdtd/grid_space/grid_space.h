#ifndef _XFDTD_CORE_GRID_SPACE_H_
#define _XFDTD_CORE_GRID_SPACE_H_

#include <xfdtd/coordinate_system/coordinate_system.h>
#include <xfdtd/exception/exception.h>
#include <xfdtd/shape/cube.h>
#include <xfdtd/shape/shape.h>

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

  double basedDx() const;

  double basedDy() const;

  double basedDz() const;

  double minDx() const;

  double minDy() const;

  double minDz() const;

  auto& gridWithMaterial() const { return _grid_with_material; }

  std::shared_ptr<GridSpace> globalGridSpace() const;

  const xt::xarray<double>& eNodeX() const;

  const xt::xarray<double>& eNodeY() const;

  const xt::xarray<double>& eNodeZ() const;

  const xt::xarray<double>& hNodeX() const;

  const xt::xarray<double>& hNodeY() const;

  const xt::xarray<double>& hNodeZ() const;

  const xt::xarray<double>& eSizeX() const;

  const xt::xarray<double>& eSizeY() const;

  const xt::xarray<double>& eSizeZ() const;

  const xt::xarray<double>& hSizeX() const;

  const xt::xarray<double>& hSizeY() const;

  const xt::xarray<double>& hSizeZ() const;

  const xt::xarray<double>& eSize(Axis::XYZ xyz) const;

  const xt::xarray<double>& hSize(Axis::XYZ xyz) const;

  std::size_t sizeX() const;

  std::size_t sizeY() const;

  std::size_t sizeZ() const;

  Grid getGrid(double x, double y, double z) const;

  Grid getGrid(const Vector& vector) const;

  Vector getGridOriginVector(const Grid& grid) const;

  Vector getGridCenterVector(const Grid& grid) const;

  Vector getGridEndVector(const Grid& grid) const;

  auto getShapeMask(const Shape* shape) const;

  auto getGridWithMaterialView(const xt::xarray<bool>& mask);

  GridBox getGridBox(const Shape* shape) const;

  virtual void correctGridSpace() = 0;

  void extendGridSpace(Axis::Direction direction, std::size_t num, double dl);

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
  GridSpace(double dx, double dy, double dz, Dimension dimension,
            xt::xarray<double> e_node_x, xt::xarray<double> e_node_y,
            xt::xarray<double> e_node_z);

  // for subGridSpace
  GridSpace(Dimension dimension, Type type, GridBox global_box, double based_dx,
            double based_dy, double based_dz, double min_dx, double min_dy,
            double min_dz, xt::xarray<double> e_node_x,
            xt::xarray<double> e_node_y, xt::xarray<double> e_node_z,
            xt::xarray<double> h_node_x, xt::xarray<double> h_node_y,
            xt::xarray<double> h_node_z, xt::xarray<double> e_size_x,
            xt::xarray<double> e_size_y, xt::xarray<double> e_size_z,
            xt::xarray<double> h_size_x, xt::xarray<double> h_size_y,
            xt::xarray<double> h_size_z);

  xt::xarray<double>& eNodeX();

  xt::xarray<double>& eNodeY();

  xt::xarray<double>& eNodeZ();

  xt::xarray<double>& hNodeX();

  xt::xarray<double>& hNodeY();

  xt::xarray<double>& hNodeZ();

  xt::xarray<double>& eSizeX();

  xt::xarray<double>& eSizeY();

  xt::xarray<double>& eSizeZ();

  xt::xarray<double>& hSizeX();

  xt::xarray<double>& hSizeY();

  xt::xarray<double>& hSizeZ();

  xt::xarray<double>& eSize(Axis::XYZ xyz);

  xt::xarray<double>& hSize(Axis::XYZ xyz);

  void setMinDx(double min_dx);

  void setMinDy(double min_dy);

  void setMinDz(double min_dz);

  virtual std::size_t handleTransformX(double x) const;

  virtual std::size_t handleTransformY(double y) const;

  virtual std::size_t handleTransformZ(double z) const;

  virtual std::size_t handleTransformXWithoutCheck(double x) const;

  virtual std::size_t handleTransformYWithoutCheck(double y) const;

  virtual std::size_t handleTransformZWithoutCheck(double z) const;

  void correctGridSpaceForOne(double dl, const xt::xarray<double>& e_node,
                              xt::xarray<double>& h_node,
                              xt::xarray<double>& e_size,
                              xt::xarray<double>& h_size);

  double _min_x, _min_y, _min_z;
  double _max_x, _max_y, _max_z;

 private:
  Dimension _dimension{Dimension::UNDEFINED};
  Type _type{Type::UNDEFINED};
  double _based_dx, _based_dy, _based_dz;
  double _min_dx, _min_dy, _min_dz;

  xt::xarray<double> _e_node_x, _e_node_y, _e_node_z;
  xt::xarray<double> _h_node_x, _h_node_y, _h_node_z;
  xt::xarray<double> _e_size_x, _e_size_y, _e_size_z;
  xt::xarray<double> _h_size_x, _h_size_y, _h_size_z;

  xt::xarray<std::shared_ptr<Grid>> _grid_with_material;

  std::weak_ptr<GridSpace> _global_grid_space;

  GridBox _global_box;

  xt::xarray<double> calculateHNode(const xt::xarray<double>& e_node);

  xt::xarray<double> calculateESize(const xt::xarray<double>& e_node);

  xt::xarray<double> calculateHSize(const xt::xarray<double>& h_node,
                                    const xt::xarray<double>& e_node,
                                    double dl);
};

inline auto GridSpace::getShapeMask(const Shape* shape) const {
  xt::xarray<bool> mask{xt::make_lambda_xfunction(
      [shape, this](auto&& g) {
        return shape->isInside(getGridCenterVector(*g));
      },
      _grid_with_material)};
  return mask;
}

inline auto GridSpace::getGridWithMaterialView(const xt::xarray<bool>& mask) {
  return xt::filter(_grid_with_material, mask);
}

}  // namespace xfdtd

#endif  // _XFDTD_CORE_GRID_SPACE_H_
