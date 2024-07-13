#ifndef __XFDTD_CORE_NFFFT_FREQUENCY_DOMAIN_H__
#define __XFDTD_CORE_NFFFT_FREQUENCY_DOMAIN_H__

#include <xfdtd/coordinate_system/coordinate_system.h>
#include <xfdtd/nffft/nffft.h>

namespace xfdtd {

class FDPlaneData;

class XFDTDNFFFTFrequencyDomainException : public XFDTDNFFFTException {
 public:
  explicit XFDTDNFFFTFrequencyDomainException(const std::string& message)
      : XFDTDNFFFTException(message) {}
};

class NFFFTFrequencyDomain : public NFFFT {
 public:
  NFFFTFrequencyDomain(Index distance_x, Index distance_y, Index distance_z,
                       Array1D<Real> frequencies);

  NFFFTFrequencyDomain(Index distance_x, Index distance_y, Index distance_z,
                       Array1D<Real> frequencies, std::string_view output_dir);

  ~NFFFTFrequencyDomain() override;

  void init(std::shared_ptr<const GridSpace> grid_space,
            std::shared_ptr<const CalculationParam> calculation_param,
            std::shared_ptr<const EMF> emf) override;

  auto initTimeDependentVariable() -> void override;

  auto update() -> void override;

  auto processFarField(const Array1D<Real>& theta, Real phi,
                       const std::string& sub_dir,
                       const Vector& origin = Vector{0.0, 0.0,
                                                     0.0}) const -> void;

  auto processFarField(Real theta, const Array1D<Real>& phi,
                       const std::string& sub_dir,
                       const Vector& origin = Vector{0.0, 0.0,
                                                     0.0}) const -> void;

  auto outputRadiationPower() -> void;

  template <Axis::Direction D, EMF::Attribute A, Axis::XYZ XYZ>
  auto equivalentSurfaceCurrent(Index freq_index) const
      -> const Array3D<std::complex<Real>>&;

  auto transformE(Index freq_index) const -> const Array1D<std::complex<Real>>&;

  auto transformH(Index freq_index) const -> const Array1D<std::complex<Real>>&;

  auto freqCount() const { return _frequencies.size(); }

 private:
  Array1D<Real> _frequencies;

  std::vector<FDPlaneData> _fd_plane_data;

  auto processFarField(const Array1D<Real>& theta, const Array1D<Real>& phi,
                       const std::string& sub_dir,
                       const Vector& origin) const -> void;

  auto generateSurface() -> void;
};

}  // namespace xfdtd

#endif  // __XFDTD_CORE_NFFFT_FREQUENCY_DOMAIN_H__
