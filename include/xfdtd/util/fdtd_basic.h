#ifndef __XFDTD_BASIC_FDTD_BASIC_H__
#define __XFDTD_BASIC_FDTD_BASIC_H__

#include <xfdtd/common/index_task.h>

#include <type_traits>

#include "xfdtd/coordinate_system/coordinate_system.h"
#include "xfdtd/electromagnetic_field/electromagnetic_field.h"

namespace xfdtd::basic {

template <typename T>
concept Number = std::is_arithmetic_v<T>;

// Index means the index of the grid space.
class GridStructure {
 public:
  static auto exSizeX(const Number auto& nx) { return nx; }

  static auto exSizeY(const Number auto& ny) { return ny + 1; }

  static auto exSizeZ(const Number auto& nz) { return nz + 1; }

  static auto exSize(const Number auto& nx, const Number auto& ny,
                     const Number auto& nz) {
    return std::array{exSizeX(nx), exSizeY(ny), exSizeZ(nz)};
  }

  static auto exFDTDUpdateXStart(const Number auto& start_x) { return start_x; }

  static auto exFDTDUpdateYStart(const Number auto& start_y) {
    return start_y + 1;
  }

  static auto exFDTDUpdateZStart(const Number auto& start_z) {
    return start_z + 1;
  }

  static auto exFDTDUpdateXEnd(const Number auto& end_x) { return end_x; }

  static auto exFDTDUpdateYEnd(const Number auto& end_y) { return end_y; }

  static auto exFDTDUpdateZEnd(const Number auto& end_z) { return end_z; }

  template <Number T>
  static auto exFDTDUpdateXRange(const Range<T>& range) {
    return makeRange(exFDTDUpdateXStart(range.start()),
                              exFDTDUpdateXEnd(range.end()));
  }

  template <Number T>
  static auto exFDTDUpdateYRange(const Range<T>& range) {
    return makeRange(exFDTDUpdateYStart(range.start()),
                              exFDTDUpdateYEnd(range.end()));
  }

  template <Number T>
  static auto exFDTDUpdateZRange(const Range<T>& range) {
    return makeRange(exFDTDUpdateZStart(range.start()),
                              exFDTDUpdateZEnd(range.end()));
  }

  template <Number T>
  static auto exFDTDUpdateTask(const Task<T>& task) {
    return makeTask(exFDTDUpdateXRange(task.xRange()),
                             exFDTDUpdateYRange(task.yRange()),
                             exFDTDUpdateZRange(task.zRange()));
  }

  static auto eySizeX(const Number auto& nx) { return nx + 1; }

  static auto eySizeY(const Number auto& ny) { return ny; }

  static auto eySizeZ(const Number auto& nz) { return nz + 1; }

  static auto eySize(const Number auto& nx, const Number auto& ny,
                     const Number auto& nz) {
    return std::array{eySizeX(nx), eySizeY(ny), eySizeZ(nz)};
  }

  static auto eyFDTDUpdateXStart(const Number auto& start_x) {
    return start_x + 1;
  }

  static auto eyFDTDUpdateYStart(const Number auto& start_y) { return start_y; }

  static auto eyFDTDUpdateZStart(const Number auto& start_z) {
    return start_z + 1;
  }

  static auto eyFDTDUpdateXEnd(const Number auto& end_x) { return end_x; }

  static auto eyFDTDUpdateYEnd(const Number auto& end_y) { return end_y; }

  static auto eyFDTDUpdateZEnd(const Number auto& end_z) { return end_z; }

  template <Number T>
  static auto eyFDTDUpdateXRange(const Range<T>& range) {
    return makeRange(eyFDTDUpdateXStart(range.start()),
                              eyFDTDUpdateXEnd(range.end()));
  }

  template <Number T>
  static auto eyFDTDUpdateYRange(const Range<T>& range) {
    return makeRange(eyFDTDUpdateYStart(range.start()),
                              eyFDTDUpdateYEnd(range.end()));
  }

  template <Number T>
  static auto eyFDTDUpdateZRange(const Range<T>& range) {
    return makeRange(eyFDTDUpdateZStart(range.start()),
                              eyFDTDUpdateZEnd(range.end()));
  }

  template <Number T>
  static auto eyFDTDUpdateTask(const Task<T>& task) {
    return makeTask(eyFDTDUpdateXRange(task.xRange()),
                             eyFDTDUpdateYRange(task.yRange()),
                             eyFDTDUpdateZRange(task.zRange()));
  }

  static auto ezSizeX(const Number auto& nx) { return nx + 1; }

  static auto ezSizeY(const Number auto& ny) { return ny + 1; }

  static auto ezSizeZ(const Number auto& nz) { return nz; }

  static auto ezSize(const Number auto& nx, const Number auto& ny,
                     const Number auto& nz) {
    return std::array{ezSizeX(nx), ezSizeY(ny), ezSizeZ(nz)};
  }

  static auto ezFDTDUpdateXStart(const Number auto& start_x) {
    return start_x + 1;
  }

  static auto ezFDTDUpdateYStart(const Number auto& start_y) {
    return start_y + 1;
  }

  static auto ezFDTDUpdateZStart(const Number auto& start_z) { return start_z; }

  static auto ezFDTDUpdateXEnd(const Number auto& end_x) { return end_x; }

  static auto ezFDTDUpdateYEnd(const Number auto& end_y) { return end_y; }

  static auto ezFDTDUpdateZEnd(const Number auto& end_z) { return end_z; }

  template <Number T>
  static auto ezFDTDUpdateXRange(const Range<T>& range) {
    return makeRange(ezFDTDUpdateXStart(range.start()),
                              ezFDTDUpdateXEnd(range.end()));
  }

  template <Number T>
  static auto ezFDTDUpdateYRange(const Range<T>& range) {
    return makeRange(ezFDTDUpdateYStart(range.start()),
                              ezFDTDUpdateYEnd(range.end()));
  }

  template <Number T>
  static auto ezFDTDUpdateZRange(const Range<T>& range) {
    return makeRange(ezFDTDUpdateZStart(range.start()),
                              ezFDTDUpdateZEnd(range.end()));
  }

  template <Number T>
  static auto ezFDTDUpdateTask(const Task<T>& task) {
    return makeTask(ezFDTDUpdateXRange(task.xRange()),
                             ezFDTDUpdateYRange(task.yRange()),
                             ezFDTDUpdateZRange(task.zRange()));
  }

  /*
  IMPORTANT:
  The H field is updated in the different order from the E field.
  Hx task is [0, nx+1), [0, ny), [0, nz)
  Hy task is [0, nx), [0, ny+1), [0, nz)
  Hz task is [0, nx), [0, ny), [0, nz+1)
  but they are updated in the same order.
  for (int k = 0; k < nz; ++k) {
    for (int j = 0; j < ny; ++j) {
      for (int i = 0; i < nx; ++i) {
        Hx...
        Hy...
        Hz...
      }
    }
  }
  The reason is that we don't care about the boundary of the grid space, so we
  don't update the Hx(nx+1, :, :) and so on even if the it can be updated.
  If hx, hy, hz can be written in the same update order, the code is simply. So
  I think it is better to write the update order in the same order.
  */
  template <Number T>
  static auto hFDTDUpdateXStart(const T& start_x) {
    return start_x;
  }

  template <Number T>
  static auto hFDTDUpdateYStart(const T& start_y) {
    return start_y;
  }

  template <Number T>
  static auto hFDTDUpdateZStart(const T& start_z) {
    return start_z;
  }

  template <Number T>
  static auto hFDTDUpdateXEnd(const T& end_x) {
    return end_x;
  }

  template <Number T>
  static auto hFDTDUpdateYEnd(const T& end_y) {
    return end_y;
  }

  template <Number T>
  static auto hFDTDUpdateZEnd(const T& end_z) {
    return end_z;
  }

  template <Number T>
  static auto hFDTDUpdateXRange(const Range<T>& range) {
    return makeRange(hFDTDUpdateXStart(range.start()),
                              hFDTDUpdateXEnd(range.end()));
  }

  template <Number T>
  static auto hFDTDUpdateYRange(const Range<T>& range) {
    return makeRange(hFDTDUpdateYStart(range.start()),
                              hFDTDUpdateYEnd(range.end()));
  }

  template <Number T>
  static auto hFDTDUpdateZRange(const Range<T>& range) {
    return makeRange(hFDTDUpdateZStart(range.start()),
                              hFDTDUpdateZEnd(range.end()));
  }

  template <Number T>
  static auto hFDTDUpdateTask(const Task<T>& task) {
    return makeTask(hFDTDUpdateXRange(task.xRange()),
                             hFDTDUpdateYRange(task.yRange()),
                             hFDTDUpdateZRange(task.zRange()));
  }

  template <Number T>
  static auto hxSizeX(const T& nx) {
    return nx + 1;
  }

  template <Number T>
  static auto hxSizeY(const T& ny) {
    return ny;
  }

  template <Number T>
  static auto hxSizeZ(const T& nz) {
    return nz;
  }

  template <Number T>
  static auto hxSize(const T& nx, const T& ny, const T& nz) {
    return std::array{hxSizeX(nx), hxSizeY(ny), hxSizeZ(nz)};
  }

  template <Number T>
  static auto hySizeX(const T& nx) {
    return nx;
  }

  template <Number T>
  static auto hySizeY(const T& ny) {
    return ny + 1;
  }

  template <Number T>
  static auto hySizeZ(const T& nz) {
    return nz;
  }

  template <Number T>
  static auto hySize(const T& nx, const T& ny, const T& nz) {
    return std::array{hySizeX(nx), hySizeY(ny), hySizeZ(nz)};
  }

  template <Number T>
  static auto hzSizeX(const T& nx) {
    return nx;
  }

  template <Number T>
  static auto hzSizeY(const T& ny) {
    return ny;
  }

  template <Number T>
  static auto hzSizeZ(const T& nz) {
    return nz + 1;
  }

  template <Number T>
  static auto hzSize(const T& nx, const T& ny, const T& nz) {
    return std::array{hzSizeX(nx), hzSizeY(ny), hzSizeZ(nz)};
  }
};

class FieldStruture {};

namespace grid_layout {

template <Number T, EMF::Attribute attribute_c, Axis::XYZ xyz_c>
inline auto diffPreNodeAxisA(T& i, T& j, T& k) {
  if constexpr (attribute_c == EMF::Attribute::E) {
    throw std::invalid_argument("not implemented");
  } else if constexpr (attribute_c == EMF::Attribute::H) {
    return;
  }
}

template <Number T, EMF::Attribute attribute_c, Axis::XYZ xyz_c>
inline auto diffPreNodeAxisB(T& i, T& j, T& k) {
  if constexpr (attribute_c == EMF::Attribute::E) {
    throw std::invalid_argument("not implemented");
  } else if constexpr (attribute_c == EMF::Attribute::H) {
    return;
  }
}

template <Number T, EMF::Attribute attribute_c, Axis::XYZ xyz_c>
inline auto diffNexNodeAxisA(T& i, T& j, T& k) {
  if constexpr (attribute_c == EMF::Attribute::E) {
    return;
  } else if constexpr (attribute_c == EMF::Attribute::H) {
    if constexpr (xyz_c == Axis::XYZ::X) {
      ++k;
    } else if constexpr (xyz_c == Axis::XYZ::Y) {
      ++i;
    } else if constexpr (xyz_c == Axis::XYZ::Z) {
      ++j;
    }
  }
}

template <Number T, EMF::Attribute attribute_c, Axis::XYZ xyz_c>
inline auto diffNexNodeAxisB(T& i, T& j, T& k) {
  if constexpr (attribute_c == EMF::Attribute::E) {
    return;
  } else if constexpr (attribute_c == EMF::Attribute::H) {
    if constexpr (xyz_c == Axis::XYZ::X) {
      ++j;
    } else if constexpr (xyz_c == Axis::XYZ::Y) {
      ++k;
    } else if constexpr (xyz_c == Axis::XYZ::Z) {
      ++i;
    }
  }
}

}  // namespace grid_layout

}  // namespace xfdtd::basic

#endif  // __XFDTD_BASIC_FDTD_BASIC_H__
