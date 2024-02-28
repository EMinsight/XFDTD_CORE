#ifndef _XFDTD_LIB_PML_H_
#define _XFDTD_LIB_PML_H_

#include <xfdtd/boundary/boundary.h>

#include <memory>
#include <xtensor/xview.hpp>

#include "corrector/corrector.h"
#include "xfdtd/coordinate_system/coordinate_system.h"

namespace xfdtd {

class XFDTDPMLException : public XFDTDBoundaryException {
 public:
  explicit XFDTDPMLException(const std::string& message = "XFDTD PML Exception")
      : XFDTDBoundaryException(message) {}
};

class PML : public Boundary {
 public:
  PML(int thickness, Axis::Direction direction, int order = 4,
      double sigma_ratio = 1, double alpha_min = 0, double alpha_max = 0.05,
      double kappa_max = 10);

  void init(double dl, double dt, std::size_t start_index, int na, int nb,
            xt::xarray<double>& ceahb, xt::xarray<double>& cebha,
            xt::xarray<double>& chaeb, xt::xarray<double>& chbea);

  PML(const PML&) = delete;

  PML(PML&&) noexcept = default;

  PML& operator=(const PML&) = delete;

  PML& operator=(PML&&) noexcept = default;

  ~PML() override = default;

  void init(std::shared_ptr<const GridSpace> grid_space,
            std::shared_ptr<CalculationParam> calculation_param,
            std::shared_ptr<EMF> emf) override;

  void correctMaterialSpace() override;

  void correctUpdateCoefficient() override;

  void correctE() override;

  void correctH() override;

  int thickness() const;

  Axis::Direction direction() const;

  Axis::XYZ subAxisA() const;

  Axis::XYZ subAxisB() const;

  Axis::XYZ mainAxis() const;

  std::size_t eNodeStartIndexMainAxis() const;

  std::size_t hNodeStartIndexMainAxis() const;

  std::size_t n() const;

  const xt::xarray<double>& eSize() const;

  const xt::xarray<double>& hSize() const;

  std::unique_ptr<Corrector> generateDomainCorrector(
      const Divider::Task<std::size_t>& task) override;

 private:
  int _thickness;
  std::size_t _n;
  Axis::Direction _direction;
  Axis::XYZ _main_axis;
  int _order;
  double _sigma_ratio;
  double _alpha_min;
  double _alpha_max;
  double _kappa_max;
  std::size_t _e_start_index, _h_start_index;
  std::size_t _na, _nb;

  xt::xarray<double> _h_size;
  xt::xarray<double> _e_size;
  xt::xarray<double> _kappa_e;
  xt::xarray<double> _kappa_h;

  xt::xarray<double> _coeff_a_e;
  xt::xarray<double> _coeff_b_e;
  xt::xarray<double> _coeff_a_h;
  xt::xarray<double> _coeff_b_h;

  xt::xarray<double> _c_ea_psi_hb;
  xt::xarray<double> _ea_psi_hb;
  xt::xarray<double> _c_eb_psi_ha;
  xt::xarray<double> _eb_psi_ha;

  xt::xarray<double> _c_ha_psi_eb;
  xt::xarray<double> _ha_psi_eb;
  xt::xarray<double> _c_hb_psi_ea;
  xt::xarray<double> _hb_psi_ea;

  bool taskContainPML(const Divider::Task<std::size_t>& task) const;

  xt::xarray<double>& eaF();

  xt::xarray<double>& ebF();

  xt::xarray<double>& haF();

  xt::xarray<double>& hbF();

  void correctCoefficientX();

  void correctCoefficientY();

  void correctCoefficientZ();

  void calRecursiveConvolutionCoeff();

  double calculateSigmaMax(double dl) const;

  xt::xarray<double> calculateRhoE(std::size_t n,
                                   const xt::xarray<double>& size) const;

  xt::xarray<double> calculateRhoM(std::size_t n,
                                   const xt::xarray<double>& size) const;

  xt::xarray<double> calculateSigma(double sigma_max,
                                    const xt::xarray<double>& rho,
                                    std::size_t order) const;

  xt::xarray<double> calculateKappa(double kappa_max,
                                    const xt::xarray<double>& rho,
                                    std::size_t order) const;

  xt::xarray<double> calculateAlpha(double alpha_min, double alpha_max,
                                    const xt::xarray<double>& rho) const;

  xt::xarray<double> calculateCoefficientA(const xt::xarray<double>& b,
                                           const xt::xarray<double>& sigma,
                                           const xt::xarray<double>& kappa,
                                           const xt::xarray<double>& alpha,
                                           const xt::xarray<double>& dl) const;

  xt::xarray<double> calculateCoefficientB(const xt::xarray<double>& sigma,
                                           const xt::xarray<double>& kappa,
                                           const xt::xarray<double>& alpha,
                                           double dt, double constant) const;

  xt::xarray<double> calculateCoeffPsi(const xt::xarray<double>& coeff,
                                       const xt::xarray<double>& kappa,
                                       const xt::xarray<double>& dl) const;
};

}  // namespace xfdtd

#endif  // _XFDTD_LIB_PML_H_
