#ifndef __XFDTD_CORE_NFFFT_FD_DATA_H__
#define __XFDTD_CORE_NFFFT_FD_DATA_H__

#include <xfdtd/calculation_param/calculation_param.h>
#include <xfdtd/coordinate_system/coordinate_system.h>
#include <xfdtd/electromagnetic_field/electromagnetic_field.h>
#include <xfdtd/grid_space/grid_space.h>
#include <xfdtd/nffft/nffft.h>
#include <xfdtd/util/transform.h>

namespace xfdtd {

class FDPlaneData {
 public:
  enum class Potential { A, F };

 public:
  FDPlaneData(std::shared_ptr<const GridSpace> grid_space,
              std::shared_ptr<const EMF> emf, Real freq,
              const IndexTask& task_xn, const IndexTask& task_xp,
              const IndexTask& task_yn, const IndexTask& task_yp,
              const IndexTask& task_zn, const IndexTask& task_zp);

  auto frequency() const -> Real;

  auto update(std::size_t current_time_step) -> void;

  auto aTheta(const Array1D<Real>& theta, const Array1D<Real>& phi,
              const Vector& origin) const -> Array1D<std::complex<Real>>;
  auto fPhi(const Array1D<Real>& theta, const Array1D<Real>& phi,
            const Vector& origin) const -> Array1D<std::complex<Real>>;

  auto aPhi(const Array1D<Real>& theta, const Array1D<Real>& phi,
            const Vector& origin) const -> Array1D<std::complex<Real>>;
  auto fTheta(const Array1D<Real>& theta, const Array1D<Real>& phi,
              const Vector& origin) const -> Array1D<std::complex<Real>>;

  auto power() const -> Real;

  auto initDFT(std::size_t total_time_step, Real dt) -> void;

  template <Potential p, transform::SCS scs>
  auto getPotential(const Array1D<Real>& theta, const Array1D<Real>& phi,
                    const Vector& origin) const -> Array1D<std::complex<Real>>;

 private:
  template <Potential potential, Axis::Direction direction>
  auto calculatePotential(const Array1D<Real>& theta, const Array1D<Real>& phi,
                          const Vector& origin,
                          const Array1D<std::complex<Real>>& transform_a,
                          const Array1D<std::complex<Real>>& transform_b) const
      -> Array1D<std::complex<Real>>;

  template <Axis::Direction direction>
  auto calculatePower() const -> std::complex<Real>;

  template <Axis::XYZ xyz>
  auto rVector(size_t i, std::size_t j, std::size_t k) const -> Vector;

  template <Axis::XYZ xyz>
  auto ds(size_t i, std::size_t j, std::size_t k) const -> Real;

 private:
  std::shared_ptr<const GridSpace> _grid_space;
  std::shared_ptr<const EMF> _emf;
  Real _freq;
  IndexTask _task_xn, _task_xp, _task_yn, _task_yp, _task_zn, _task_zp;
  Array3D<std::complex<Real>> _jx_yn, _jx_yp, _jx_zn, _jx_zp;
  Array3D<std::complex<Real>> _jy_xn, _jy_xp, _jy_zn, _jy_zp;
  Array3D<std::complex<Real>> _jz_xn, _jz_xp, _jz_yn, _jz_yp;
  Array3D<std::complex<Real>> _mx_yn, _mx_yp, _mx_zn, _mx_zp;
  Array3D<std::complex<Real>> _my_xn, _my_xp, _my_zn, _my_zp;
  Array3D<std::complex<Real>> _mz_xn, _mz_xp, _mz_yn, _mz_yp;

  Array1D<std::complex<Real>> _transform_e;
  Array1D<std::complex<Real>> _transform_h;

  template <Axis::Direction direction>
  auto calculateJ(std::size_t current_time_step) -> void;

  template <Axis::Direction direction>
  auto calculateM(std::size_t current_time_step) -> void;

  template <Axis::Direction direction>
  auto task() const -> const IndexTask&;

  template <Axis::XYZ xyz>
  static constexpr auto tangentialHaFieldEnum() -> EMF::Field;

  template <Axis::XYZ xyz>
  static constexpr auto tangentialHbFieldEnum() -> EMF::Field;

  template <Axis::XYZ xyz>
  static constexpr auto tangentialEaFieldEnum() -> EMF::Field;

  template <Axis::XYZ xyz>
  static constexpr auto tangentialEbFieldEnum() -> EMF::Field;

  template <Axis::XYZ xyz, EMF::Field f>
  auto interpolate(const std::size_t& i, const std::size_t& j,
                   const std::size_t& k) -> Real;

  template <Potential potential, Axis::Direction direction>
  auto surfaceCurrent()
      -> std::tuple<Array3D<std::complex<Real>>&, Array3D<std::complex<Real>>&>;

  template <Axis::Direction direction>
  auto surfaceJ()
      -> std::tuple<Array3D<std::complex<Real>>&, Array3D<std::complex<Real>>&>;

  template <Axis::Direction direction>
  auto surfaceM()
      -> std::tuple<Array3D<std::complex<Real>>&, Array3D<std::complex<Real>>&>;

  template <Potential potential, Axis::Direction direction>
  auto surfaceCurrent() const -> std::tuple<const Array3D<std::complex<Real>>&,
                                            const Array3D<std::complex<Real>>&>;

  template <Axis::Direction direction>
  auto surfaceJ() const -> std::tuple<const Array3D<std::complex<Real>>&,
                                      const Array3D<std::complex<Real>>&>;

  template <Axis::Direction direction>
  auto surfaceM() const -> std::tuple<const Array3D<std::complex<Real>>&,
                                      const Array3D<std::complex<Real>>&>;
};

template <Axis::Direction direction>
inline auto FDPlaneData::task() const -> const IndexTask& {
  if constexpr (direction == Axis::Direction::XN) {
    return _task_xn;
  } else if constexpr (direction == Axis::Direction::XP) {
    return _task_xp;
  } else if constexpr (direction == Axis::Direction::YN) {
    return _task_yn;
  } else if constexpr (direction == Axis::Direction::YP) {
    return _task_yp;
  } else if constexpr (direction == Axis::Direction::ZN) {
    return _task_zn;
  } else if constexpr (direction == Axis::Direction::ZP) {
    return _task_zp;
  } else {
    throw XFDTDNFFFTException("Invalid direction");
  }
}

template <FDPlaneData::Potential potential, Axis::Direction direction>
inline auto FDPlaneData::surfaceCurrent()
    -> std::tuple<Array3D<std::complex<Real>>&, Array3D<std::complex<Real>>&> {
  if constexpr (potential == Potential::A) {
    return surfaceJ<direction>();
  } else if constexpr (potential == Potential::F) {
    return surfaceM<direction>();
  } else {
    throw XFDTDNFFFTException("Invalid potential");
  }
}

template <Axis::Direction direction>
inline auto FDPlaneData::surfaceJ()
    -> std::tuple<Array3D<std::complex<Real>>&, Array3D<std::complex<Real>>&> {
  if constexpr (direction == Axis::Direction::XN) {
    return std::make_tuple(std::ref(_jy_xn), std::ref(_jz_xn));
  } else if constexpr (direction == Axis::Direction::XP) {
    return std::make_tuple(std::ref(_jy_xp), std::ref(_jz_xp));
  } else if constexpr (direction == Axis::Direction::YN) {
    return std::make_tuple(std::ref(_jz_yn), std::ref(_jx_yn));
  } else if constexpr (direction == Axis::Direction::YP) {
    return std::make_tuple(std::ref(_jz_yp), std::ref(_jx_yp));
  } else if constexpr (direction == Axis::Direction::ZN) {
    return std::make_tuple(std::ref(_jx_zn), std::ref(_jy_zn));
  } else if constexpr (direction == Axis::Direction::ZP) {
    return std::make_tuple(std::ref(_jx_zp), std::ref(_jy_zp));
  } else {
    throw XFDTDNFFFTException("Invalid direction");
  }
}

template <Axis::Direction direction>
inline auto FDPlaneData::surfaceM()
    -> std::tuple<Array3D<std::complex<Real>>&, Array3D<std::complex<Real>>&> {
  if constexpr (direction == Axis::Direction::XN) {
    return std::make_tuple(std::ref(_my_xn), std::ref(_mz_xn));
  } else if constexpr (direction == Axis::Direction::XP) {
    return std::make_tuple(std::ref(_my_xp), std::ref(_mz_xp));
  } else if constexpr (direction == Axis::Direction::YN) {
    return std::make_tuple(std::ref(_mz_yn), std::ref(_mx_yn));
  } else if constexpr (direction == Axis::Direction::YP) {
    return std::make_tuple(std::ref(_mz_yp), std::ref(_mx_yp));
  } else if constexpr (direction == Axis::Direction::ZN) {
    return std::make_tuple(std::ref(_mx_zn), std::ref(_my_zn));
  } else if constexpr (direction == Axis::Direction::ZP) {
    return std::make_tuple(std::ref(_mx_zp), std::ref(_my_zp));
  } else {
    throw XFDTDNFFFTException("Invalid direction");
  }
}

template <FDPlaneData::Potential potential, Axis::Direction direction>
inline auto FDPlaneData::surfaceCurrent() const
    -> std::tuple<const Array3D<std::complex<Real>>&,
                  const Array3D<std::complex<Real>>&> {
  if constexpr (potential == Potential::A) {
    return surfaceJ<direction>();
  } else if constexpr (potential == Potential::F) {
    return surfaceM<direction>();
  } else {
    throw XFDTDNFFFTException("Invalid potential");
  }
}

template <Axis::Direction direction>
inline auto FDPlaneData::surfaceJ() const
    -> std::tuple<const Array3D<std::complex<Real>>&,
                  const Array3D<std::complex<Real>>&> {
  if constexpr (direction == Axis::Direction::XN) {
    return std::make_tuple(std::ref(_jy_xn), std::ref(_jz_xn));
  } else if constexpr (direction == Axis::Direction::XP) {
    return std::make_tuple(std::ref(_jy_xp), std::ref(_jz_xp));
  } else if constexpr (direction == Axis::Direction::YN) {
    return std::make_tuple(std::ref(_jz_yn), std::ref(_jx_yn));
  } else if constexpr (direction == Axis::Direction::YP) {
    return std::make_tuple(std::ref(_jz_yp), std::ref(_jx_yp));
  } else if constexpr (direction == Axis::Direction::ZN) {
    return std::make_tuple(std::ref(_jx_zn), std::ref(_jy_zn));
  } else if constexpr (direction == Axis::Direction::ZP) {
    return std::make_tuple(std::ref(_jx_zp), std::ref(_jy_zp));
  } else {
    throw XFDTDNFFFTException("Invalid direction");
  }
}

template <Axis::Direction direction>
inline auto FDPlaneData::surfaceM() const
    -> std::tuple<const Array3D<std::complex<Real>>&,
                  const Array3D<std::complex<Real>>&> {
  if constexpr (direction == Axis::Direction::XN) {
    return std::make_tuple(std::ref(_my_xn), std::ref(_mz_xn));
  } else if constexpr (direction == Axis::Direction::XP) {
    return std::make_tuple(std::ref(_my_xp), std::ref(_mz_xp));
  } else if constexpr (direction == Axis::Direction::YN) {
    return std::make_tuple(std::ref(_mz_yn), std::ref(_mx_yn));
  } else if constexpr (direction == Axis::Direction::YP) {
    return std::make_tuple(std::ref(_mz_yp), std::ref(_mx_yp));
  } else if constexpr (direction == Axis::Direction::ZN) {
    return std::make_tuple(std::ref(_mx_zn), std::ref(_my_zn));
  } else if constexpr (direction == Axis::Direction::ZP) {
    return std::make_tuple(std::ref(_mx_zp), std::ref(_my_zp));
  } else {
    throw XFDTDNFFFTException("Invalid direction");
  }
}

}  // namespace xfdtd

#endif  // __XFDTD_CORE_NFFFT_FD_DATA_H__
