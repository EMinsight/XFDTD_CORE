#ifndef __XFDTD_CORE_GRID_BOX_H__
#define __XFDTD_CORE_GRID_BOX_H__

#include "xfdtd/grid_space/grid.h"

namespace xfdtd {

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

}  // namespace xfdtd

#endif  // __XFDTD_CORE_GRID_BOX_H__
