#ifndef __XFDTD_CORE_GRID_SPACE_DATA_H__
#define __XFDTD_CORE_GRID_SPACE_DATA_H__

namespace xfdtd {
template <typename Arr1D, typename Index, typename Real>
class GridSpaceData {
  friend class GridSpaceHD;

 public:
  // GridSpaceData() = default;

  // GridSpaceData(Real based_dx, Real based_dy, Real based_dz, Real min_dx,
  //               Real min_dy, Real min_dz, Arr1D *e_node_x, Arr1D *e_node_y,
  //               Arr1D *e_node_z, Arr1D *h_node_x, Arr1D *h_node_y,
  //               Arr1D *h_node_z, Arr1D *e_size_x, Arr1D *e_size_y,
  //               Arr1D *e_size_z, Arr1D *h_size_x, Arr1D *h_size_y,
  //               Arr1D *h_size_z);

  // /**
  //  * @brief e node means that the corner of the cell
  //  *
  //  * @return const Arr1D&
  //  */
  // auto eNodeX() const -> const Arr1D &;

  // /**
  //  * @brief h node means that the center of the cell
  //  *
  //  * @return const Arr1D&
  //  */
  // auto hNodeX() const -> const Arr1D &;

  // auto eNodeY() const -> const Arr1D &;

  // auto eNodeZ() const -> const Arr1D &;

  // auto hNodeY() const -> const Arr1D &;

  // auto hNodeZ() const -> const Arr1D &;

  // auto eSizeX() const -> const Arr1D &;

  // auto eSizeY() const -> const Arr1D &;

  // auto eSizeZ() const -> const Arr1D &;

  // auto hSizeX() const -> const Arr1D &;

  // auto hSizeY() const -> const Arr1D &;

  // auto hSizeZ() const -> const Arr1D &;

  // /**
  //  * @brief The number of cells in x direction
  //  *
  //  * @return Index
  //  */
  // auto sizeX() const -> Index;

  // /**
  //  * @brief The number of cells in y direction
  //  *
  //  * @return Index
  //  */
  // auto sizeY() const -> Index;

  // /**
  //  * @brief The number of cells in z direction
  //  *
  //  * @return Index
  //  */
  // auto sizeZ() const -> Index;


  Real _based_dx{}, _based_dy{}, _based_dz{};
  Real _min_dx{}, _min_dy{}, _min_dz{};

  const Arr1D *_e_node_x{}, *_e_node_y{}, *_e_node_z{};
  const Arr1D *_h_node_x{}, *_h_node_y{}, *_h_node_z{};
  const Arr1D *_e_size_x{}, *_e_size_y{}, *_e_size_z{};
  const Arr1D *_h_size_x{}, *_h_size_y{}, *_h_size_z{};
};

}  // namespace xfdtd

#endif  // __XFDTD_CORE_GRID_SPACE_DATA_H__
