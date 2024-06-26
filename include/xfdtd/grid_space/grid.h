#ifndef __XFDTD_CORE_GRID_H__
#define __XFDTD_CORE_GRID_H__

#include "xfdtd/common/type_define.h"

namespace xfdtd {
class Grid {
 public:
  Grid() = default;

  Grid(Index i, Index j, Index k, Index material_index = -1);

  Grid(const Grid&) = default;

  Grid(Grid&&) noexcept = default;

  Grid& operator=(const Grid&) = default;

  Grid& operator=(Grid&&) noexcept = default;

  ~Grid() = default;

  Grid operator+(const Grid& grid) const;

  Grid operator-(const Grid& grid) const;

  Index i() const;

  Index j() const;

  Index k() const;

  Index materialIndex() const;

  void setMaterialIndex(Index index);

  std::string toString() const;

 private:
  Index _i{0}, _j{0}, _k{0}, _material_index{static_cast<size_t>(-1)};
};
}  // namespace xfdtd

#endif  // __XFDTD_CORE_GRID_H__
