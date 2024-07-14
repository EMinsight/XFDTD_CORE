#ifndef __XFDTD_CORE_NFFFT_TD_H__
#define __XFDTD_CORE_NFFFT_TD_H__

#include <xfdtd/nffft/nffft.h>

#include "xfdtd/common/index_task.h"
#include "xfdtd/common/type_define.h"
#include "xfdtd/coordinate_system/coordinate_system.h"
#include "xfdtd/electromagnetic_field/electromagnetic_field.h"

namespace xfdtd {

class XFDTDNFFFTTimeDomainException : public XFDTDNFFFTException {
 public:
  explicit XFDTDNFFFTTimeDomainException(const std::string& message)
      : XFDTDNFFFTException(message) {}
};

template <Axis::Direction D>
class TDPlaneData;

/**
 * @brief Near Field Far Field Transform in Time Domain
 *
 */
class NFFFTTimeDomain : public NFFFT {
 public:
  NFFFTTimeDomain(Index distance_x, Index distance_y, Index distance_z,
                  Real theta, Real phi);

  NFFFTTimeDomain(Index distance_x, Index distance_y, Index distance_z,
                  std::string_view output_dir, Real theta, Real phi);

  ~NFFFTTimeDomain() override;

  auto init(std::shared_ptr<const GridSpace> grid_space,
            std::shared_ptr<const CalculationParam> calculation_param,
            std::shared_ptr<const EMF> emf) -> void override;

  auto initTimeDependentVariable() -> void override;

  auto update() -> void override;

  auto processFarField() const -> void;

  auto observationDirection() const -> Vector;

  auto wx() const -> Array1D<Real>;

  auto wy() const -> Array1D<Real>;

  auto wz() const -> Array1D<Real>;

  auto ux() const -> Array1D<Real>;

  auto uy() const -> Array1D<Real>;

  auto uz() const -> Array1D<Real>;

  template <Axis::Direction D, EMF::Attribute A, Axis::XYZ XYZ>
  auto equivalentSurfaceCurrent() const -> const Array1D<Real>&;

  template <Axis::Direction D, EMF::Attribute A, Axis::XYZ XYZ>
  auto fieldPrev() const -> const Array2D<Real>&;

  template <EMF::Attribute A>
  auto distanceRange() const -> Range<Real>;

 private:
  Real _theta, _phi;

  std::shared_ptr<TDPlaneData<Axis::Direction::XN>> _td_plane_xn;
  std::shared_ptr<TDPlaneData<Axis::Direction::XP>> _td_plane_xp;
  std::shared_ptr<TDPlaneData<Axis::Direction::YN>> _td_plane_yn;
  std::shared_ptr<TDPlaneData<Axis::Direction::YP>> _td_plane_yp;
  std::shared_ptr<TDPlaneData<Axis::Direction::ZN>> _td_plane_zn;
  std::shared_ptr<TDPlaneData<Axis::Direction::ZP>> _td_plane_zp;

  auto cosTheta() const -> Real;

  auto sinTheta() const -> Real;

  auto cosPhi() const -> Real;

  auto sinPhi() const -> Real;
};

}  // namespace xfdtd

#endif  // __XFDTD_CORE_NFFFT_TD_H__
