#ifndef __XFDTD_CORE_DISPERSIVE_MATERIAL_EQUATION_H__
#define __XFDTD_CORE_DISPERSIVE_MATERIAL_EQUATION_H__

#include <xfdtd/common/constant.h>
#include <xfdtd/common/type_define.h>
#include <xfdtd/material/dispersive_material.h>

namespace xfdtd {

class XFDTDLinearDispersiveMaterialEquationException
    : public XFDTDLinearDispersiveMaterialException {
 public:
  explicit XFDTDLinearDispersiveMaterialEquationException(
      const std::string& message)
      : XFDTDLinearDispersiveMaterialException(message){};
};

class LinearDispersiveMaterialEquation {
 public:
  virtual ~LinearDispersiveMaterialEquation() = default;

  virtual auto numPoles() const -> Index = 0;

  virtual auto susceptibility(Real freq, Index p) const
      -> std::complex<Real> = 0;
};

/**
 * @brief the relative permittivity described by:
 /begin{equation}
    \varepsilon_r(\omega) = \varepsilon_{\infty} + \sum_{p=0}^{N-1}
    \frac{a_{1p}j\omega + a_{0p}}{b_{2p}(j\omega)^2 + b_{1p}j\omega + b_{0p}}
 /end{equation}
 *
 */
class MLorentzEqDecision : public LinearDispersiveMaterialEquation {
 public:
  MLorentzEqDecision(Array1D<Real> a_0, Array1D<Real> a_1, Array1D<Real> b_0,
                     Array1D<Real> b_1, Array1D<Real> b_2);

  ~MLorentzEqDecision() override = default;

  auto a0() const { return _a_0; }

  auto a1() const { return _a_1; }

  auto b0() const { return _b_0; }

  auto b1() const { return _b_1; }

  auto b2() const { return _b_2; }

  auto numPoles() const -> Index override { return _num_p; }

  auto susceptibility(Real freq, Index p) const -> std::complex<Real> override;

 private:
  Array1D<Real> _a_0;
  Array1D<Real> _a_1;
  Array1D<Real> _b_0;
  Array1D<Real> _b_1;
  Array1D<Real> _b_2;
  Index _num_p;
};

class DebyeEqDecision : public MLorentzEqDecision {
 public:
  DebyeEqDecision(Real epsilon_inf, Array1D<Real> epsilon_static,
                  Array1D<Real> tau);

  auto epsilonStatic() const { return _epsilon_static; }

  auto deltaEpsilon() const { return a0(); }

  auto epsilonInf() const { return epsilonStatic() - deltaEpsilon(); }

  auto tau() const { return b1(); }

  ~DebyeEqDecision() override = default;

 private:
  Array1D<Real> _epsilon_static;
};

class DrudeEqDecision : public MLorentzEqDecision {
 public:
  DrudeEqDecision(Array1D<Real> omega_p, Array1D<Real> gamma);

  ~DrudeEqDecision() override = default;

  auto omegaP() const { return xt::sqrt(a0()); }

  auto gamma() const { return b1(); }
};

}  // namespace xfdtd

#endif  // __XFDTD_CORE_DISPERSIVE_MATERIAL_EQUATION_H__
