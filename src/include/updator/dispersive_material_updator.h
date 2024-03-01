#ifndef _XFDTD_LIB_DISPERSIVE_MATERIAL_UPDATOR_H_
#define _XFDTD_LIB_DISPERSIVE_MATERIAL_UPDATOR_H_

#include <cstddef>
#include <memory>
#include <unordered_map>
#include <vector>
#include <xtensor/xtensor.hpp>

#include "material/dispersive_solver_common.h"
#include "updator/basic_updator.h"
#include "xfdtd/material/dispersive_material.h"

namespace xfdtd {

// class Object;

/**
 * @brief using ADE-FDTD
 *
 */
class LinearDispersiveMaterialADEUpdator : public BasicUpdator {
 public:
  LinearDispersiveMaterialADEUpdator(
      std::shared_ptr<const GridSpace> grid_space,
      std::shared_ptr<const CalculationParam> calculation_param,
      std::shared_ptr<EMF> emf, Divider::IndexTask task);

  LinearDispersiveMaterialADEUpdator(
      const LinearDispersiveMaterialADEUpdator&) = delete;

  LinearDispersiveMaterialADEUpdator(
      LinearDispersiveMaterialADEUpdator&&) noexcept = default;

  LinearDispersiveMaterialADEUpdator& operator=(
      const LinearDispersiveMaterialADEUpdator&) = delete;

  LinearDispersiveMaterialADEUpdator& operator=(
      LinearDispersiveMaterialADEUpdator&&) noexcept = default;

  ~LinearDispersiveMaterialADEUpdator() override = default;

 protected:
  //   auto linearDispersiveMaterialPtr() const { return
  //   _dispersive_martial.get(); }

 private:
  //   std::shared_ptr<LinearDispersiveMaterial> _dispersive_martial;
};

class LorentzADEUpdator : public LinearDispersiveMaterialADEUpdator {
 public:
  LorentzADEUpdator(std::shared_ptr<const GridSpace> grid_space,
                    std::shared_ptr<const CalculationParam> calculation_param,
                    std::shared_ptr<EMF> emf, Divider::IndexTask task);

  LorentzADEUpdator(const LorentzADEUpdator&) = delete;

  LorentzADEUpdator(LorentzADEUpdator&&) noexcept = default;

  LorentzADEUpdator& operator=(const LorentzADEUpdator&) = delete;

  LorentzADEUpdator& operator=(LorentzADEUpdator&&) noexcept = delete;

  ~LorentzADEUpdator() override = default;

  void updateE() override;

 protected:
 private:
  struct LorentzADERecord {
    xt::xtensor<double, 4> _jx_prev, _jy_prev, _jz_prev, _jx, _jy, _jz;
  };
  std::unordered_map<std::size_t, std::size_t> _lorentz_map;
  std::vector<std::shared_ptr<LorentzMedium>> _lorentz_mediums;
  std::vector<ade::LorentzCoeff> _coeff;
  std::vector<LorentzADERecord> _j_record;

  void init();

  void allocate();
};

class DrudeADEUpdator : public LinearDispersiveMaterialADEUpdator {
 public:
  DrudeADEUpdator(std::shared_ptr<const GridSpace> grid_space,
                  std::shared_ptr<const CalculationParam> calculation_param,
                  std::shared_ptr<EMF> emf, Divider::IndexTask task);

  DrudeADEUpdator(const DrudeADEUpdator&) = delete;

  DrudeADEUpdator(DrudeADEUpdator&&) noexcept = default;

  DrudeADEUpdator& operator=(const DrudeADEUpdator&) = delete;

  DrudeADEUpdator& operator=(DrudeADEUpdator&&) noexcept = default;

  ~DrudeADEUpdator() override = default;

  void updateE() override;

 private:
  struct DrudeADERecord {
    xt::xtensor<double, 4> _jx, _jy, _jz;
  };

  std::unordered_map<std::size_t, std::size_t> _drude_map;
  std::vector<std::shared_ptr<DrudeMedium>> _drude_mediums;
  std::vector<ade::DrudeCoeff> _coeff;
  std::vector<DrudeADERecord> _j_record;

  void init();

  void allocate();
};

class DebyeADEUpdator : public LinearDispersiveMaterialADEUpdator {
 public:
  DebyeADEUpdator(std::shared_ptr<const GridSpace> grid_space,
                  std::shared_ptr<const CalculationParam> calculation_param,
                  std::shared_ptr<EMF> emf, Divider::IndexTask task);

  DebyeADEUpdator(const DebyeADEUpdator&) = delete;

  DebyeADEUpdator(DebyeADEUpdator&&) noexcept = default;

  DebyeADEUpdator& operator=(const DebyeADEUpdator&) = delete;

  DebyeADEUpdator& operator=(DebyeADEUpdator&&) noexcept = default;

  ~DebyeADEUpdator() override = default;

  void updateE() override;

 private:
  struct DebyeADERecord {
    xt::xtensor<double, 4> _jx, _jy, _jz;
  };

  std::unordered_map<std::size_t, std::size_t> _debye_map;
  std::vector<std::shared_ptr<DebyeMedium>> _debye_mediums;
  std::vector<ade::DebyCoeff> _coeff;
  std::vector<DebyeADERecord> _j_record;

  void init();
};

}  // namespace xfdtd

#endif  // _XFDTD_LIB_DISPERSIVE_MATERIAL_UPDATOR_H_
