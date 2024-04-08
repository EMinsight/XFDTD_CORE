#ifndef _XFDTD_CORE_PML_H_
#define _XFDTD_CORE_PML_H_

#include <xfdtd/boundary/boundary.h>
#include <xfdtd/common/index_task.h>
#include <xfdtd/common/type_define.h>
#include <xfdtd/coordinate_system/coordinate_system.h>

namespace xfdtd {

class XFDTDPMLException : public XFDTDBoundaryException {
 public:
  explicit XFDTDPMLException(const std::string& message = "XFDTD PML Exception")
      : XFDTDBoundaryException(message) {}
};

class PML : public Boundary {
 public:
  PML(int thickness, Axis::Direction direction, int order = 4,
      Real sigma_ratio = 1, Real alpha_min = 0, Real alpha_max = 0.05,
      Real kappa_max = 10);

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

  int thickness() const;

  Axis::Direction direction() const;

  Axis::XYZ subAxisA() const;

  Axis::XYZ subAxisB() const;

  Axis::XYZ mainAxis() const;

  std::size_t globalENodeStartIndexMainAxis() const;

  std::size_t globalHNodeStartIndexMainAxis() const;

  std::size_t nodeENodeStartIndexMainAxis() const;

  std::size_t nodeHNodeStartIndexMainAxis() const;

  std::size_t n() const;

  std::size_t nodeN() const;

  const Array1D<Real>& globalESize() const;

  const Array1D<Real>& globalHSize() const;

  std::unique_ptr<Corrector> generateDomainCorrector(
      const Task<std::size_t>& task) override;

 private:
  int _thickness;
  std::size_t _n;
  Axis::Direction _direction;
  Axis::XYZ _main_axis;
  int _order;
  Real _sigma_ratio;
  Real _alpha_min;
  Real _alpha_max;
  Real _kappa_max;
  std::size_t _global_e_start_index, _global_h_start_index;
  std::size_t _global_na, _global_nb;
  IndexTask _pml_global_task_abc;
  IndexTask _pml_global_task;
  IndexTask _pml_node_task_abc;
  IndexTask _pml_node_task;

  std::size_t _node_e_start_index, _node_h_start_index;

  Array1D<Real> _global_h_size;
  Array1D<Real> _global_e_size;
  Array1D<Real> _kappa_e;
  Array1D<Real> _kappa_h;

  Array1D<Real> _coeff_a_e;
  Array1D<Real> _coeff_b_e;
  Array1D<Real> _coeff_a_h;
  Array1D<Real> _coeff_b_h;

  Array3D<Real> _c_ea_psi_hb;
  Array3D<Real> _ea_psi_hb;
  Array3D<Real> _c_eb_psi_ha;
  Array3D<Real> _eb_psi_ha;

  Array3D<Real> _c_ha_psi_eb;
  Array3D<Real> _ha_psi_eb;
  Array3D<Real> _c_hb_psi_ea;
  Array3D<Real> _hb_psi_ea;

  //   bool taskContainPML(const Task<std::size_t>& task) const;

  Array3D<Real>& eaF();

  Array3D<Real>& ebF();

  Array3D<Real>& haF();

  Array3D<Real>& hbF();

  void correctCoefficientX();

  void correctCoefficientY();

  void correctCoefficientZ();

  void calRecursiveConvolutionCoeff();

  Real calculateSigmaMax(Real dl) const;

  Array1D<Real> calculateRhoE(std::size_t n,
                                   const Array1D<Real>& size) const;

  Array1D<Real> calculateRhoM(std::size_t n,
                                   const Array1D<Real>& size) const;

  Array1D<Real> calculateSigma(Real sigma_max,
                                    const Array1D<Real>& rho,
                                    std::size_t order) const;

  Array1D<Real> calculateKappa(Real kappa_max,
                                    const Array1D<Real>& rho,
                                    std::size_t order) const;

  Array1D<Real> calculateAlpha(Real alpha_min, Real alpha_max,
                                    const Array1D<Real>& rho) const;

  Array1D<Real> calculateCoefficientA(const Array1D<Real>& b,
                                           const Array1D<Real>& sigma,
                                           const Array1D<Real>& kappa,
                                           const Array1D<Real>& alpha,
                                           const Array1D<Real>& dl) const;

  Array1D<Real> calculateCoefficientB(const Array1D<Real>& sigma,
                                           const Array1D<Real>& kappa,
                                           const Array1D<Real>& alpha,
                                           Real dt, Real constant) const;

  Array1D<Real> calculateCoeffPsi(const Array1D<Real>& coeff,
                                       const Array1D<Real>& kappa,
                                       const Array1D<Real>& dl) const;
};

}  // namespace xfdtd

#endif  // _XFDTD_CORE_PML_H_
