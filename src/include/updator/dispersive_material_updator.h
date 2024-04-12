#ifndef _XFDTD_CORE_DISPERSIVE_MATERIAL_UPDATOR_H_
#define _XFDTD_CORE_DISPERSIVE_MATERIAL_UPDATOR_H_

#include <updator/basic_updator.h>
#include <xfdtd/calculation_param/calculation_param.h>
#include <xfdtd/common/index_task.h>
#include <xfdtd/common/type_define.h>
#include <xfdtd/exception/exception.h>
#include <xfdtd/material/dispersive_material.h>

#include <memory>
#include <vector>

namespace xfdtd {

class XFDTDLinearDispersiveMaterialException : public XFDTDException {
 public:
  explicit XFDTDLinearDispersiveMaterialException(const std::string& message)
      : XFDTDException(message) {}
};

/**
 * @brief using ADE-FDTD
 *
 */
class LinearDispersiveMaterialADEUpdator : public BasicUpdator3D {
 public:
  class ADECorrector {
   public:
    ADECorrector(Index num_pole, IndexTask task, Real coeff_j,
                 std::shared_ptr<const CalculationParam> calculation_param,
                 std::shared_ptr<EMF> emf);

    virtual ~ADECorrector() = default;

    virtual auto updateEx(Index node_i, Index node_j, Index node_k) -> void = 0;

    virtual auto updateEy(Index node_i, Index node_j, Index node_k) -> void = 0;

    virtual auto updateEz(Index node_i, Index node_j, Index node_k) -> void = 0;

    auto numPole() const -> Index { return _num_pole; }

    auto task() const -> IndexTask { return _task; }

    auto coeffJ() const -> Real { return _coeff_j; }

   protected:
    Array4D<Real> _jx, _jy, _jz;
    std::shared_ptr<const CalculationParam> _calculation_param;
    std::shared_ptr<EMF> _emf;

   private:
    Index _num_pole;
    IndexTask _task;
    Real _coeff_j;
  };

 public:
  LinearDispersiveMaterialADEUpdator(
      std::vector<std::shared_ptr<Material>> material_arr,
      std::shared_ptr<const GridSpace> grid_space,
      std::shared_ptr<const CalculationParam> calculation_param,
      std::shared_ptr<EMF> emf, IndexTask task);

  ~LinearDispersiveMaterialADEUpdator() override = default;

  auto updateE() -> void override;

 protected:
  auto updateEEdge() -> void override;

 private:
  std::vector<Index> _map;
  std::vector<std::unique_ptr<ADECorrector>> _ade_correctors;

  auto init(std::vector<std::shared_ptr<Material>> material_arr) -> void;
};

class DebyeADECorrector
    : public LinearDispersiveMaterialADEUpdator::ADECorrector {
 public:
  DebyeADECorrector(const DebyeMedium& debye_medium,
                    std::shared_ptr<const CalculationParam> calculation_param,
                    std::shared_ptr<EMF> emf, IndexTask task);

  ~DebyeADECorrector() override = default;

  auto updateEx(Index node_i, Index node_j, Index node_k) -> void override;

  auto updateEy(Index node_i, Index node_j, Index node_k) -> void override;

  auto updateEz(Index node_i, Index node_j, Index node_k) -> void override;

 protected:
  auto calculateJSum(Index node_i, Index node_j, Index node_k,
                     const Array4D<Real>& j_arr) const -> Real;

  auto updateJ(Index node_i, Index node_j, Index node_k, Real e_next,
               Real e_cur, Real dt, Array4D<Real>& j_arr) -> void;

 private:
  Array1D<Real> _k;
  Array1D<Real> _beta;
};

class DrudeADECorrector
    : public LinearDispersiveMaterialADEUpdator::ADECorrector {
 public:
  DrudeADECorrector(const DrudeMedium& drude_medium,
                    std::shared_ptr<const CalculationParam> calculation_param,
                    std::shared_ptr<EMF> emf, IndexTask task);

  ~DrudeADECorrector() override = default;

  auto updateEx(Index node_i, Index node_j, Index node_k) -> void override;

  auto updateEy(Index node_i, Index node_j, Index node_k) -> void override;

  auto updateEz(Index node_i, Index node_j, Index node_k) -> void override;

 protected:
  auto calculateJSum(Index node_i, Index node_j, Index node_k,
                     const Array4D<Real>& j_arr) const -> Real;

  auto updateJ(Index node_i, Index node_j, Index node_k, Real e_next,
               Real e_cur, Array4D<Real>& j_arr) -> void;

 private:
  Array1D<Real> _k;
  Array1D<Real> _beta;
};

class LorentzADECorrector
    : public LinearDispersiveMaterialADEUpdator::ADECorrector {
 public:
  LorentzADECorrector(const LorentzMedium& lorentz_medium,
                      std::shared_ptr<const CalculationParam> calculation_param,
                      std::shared_ptr<EMF> emf, IndexTask task);

  ~LorentzADECorrector() override = default;

  auto updateEx(Index node_i, Index node_j, Index node_k) -> void override;

  auto updateEy(Index node_i, Index node_j, Index node_k) -> void override;

  auto updateEz(Index node_i, Index node_j, Index node_k) -> void override;

 private:
  Array4D<Real> _jx_prev, _jy_prev, _jz_prev;
  Array3D<Real> _ex_prev, _ey_prev, _ez_prev;

  Array1D<Real> _alpha, _xi, _gamma;
  Real _c1;

  auto calculateJSum(Index node_i, Index node_j, Index node_k,
                     const Array4D<Real>& j_arr,
                     const Array4D<Real>& j_prev_arr) const -> Real;

  auto updateJ(Index node_i, Index node_j, Index node_k, Real e_next,
               Real e_prev, Real dt, Array4D<Real>& j_arr,
               Array4D<Real>& j_prev_arr) -> void;
};

}  // namespace xfdtd

#endif  // _XFDTD_CORE_DISPERSIVE_MATERIAL_UPDATOR_H_
