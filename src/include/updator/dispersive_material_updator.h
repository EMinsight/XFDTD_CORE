#ifndef _XFDTD_CORE_DISPERSIVE_MATERIAL_UPDATOR_H_
#define _XFDTD_CORE_DISPERSIVE_MATERIAL_UPDATOR_H_

#include <xfdtd/material/dispersive_material.h>

#include <memory>
#include <unordered_map>
#include <vector>

#include "updator/basic_updator.h"

namespace xfdtd {

// class Object;

/**
 * @brief using ADE-FDTD
 *
 */
class LinearDispersiveMaterialADEUpdator : public BasicUpdator3D {
 public:
  LinearDispersiveMaterialADEUpdator(
      std::shared_ptr<const GridSpace> grid_space,
      std::shared_ptr<const CalculationParam> calculation_param,
      std::shared_ptr<EMF> emf, IndexTask task);

  LinearDispersiveMaterialADEUpdator(
      const LinearDispersiveMaterialADEUpdator&) = delete;

  LinearDispersiveMaterialADEUpdator(
      LinearDispersiveMaterialADEUpdator&&) noexcept = default;

  LinearDispersiveMaterialADEUpdator& operator=(
      const LinearDispersiveMaterialADEUpdator&) = delete;

  LinearDispersiveMaterialADEUpdator& operator=(
      LinearDispersiveMaterialADEUpdator&&) noexcept = default;

  ~LinearDispersiveMaterialADEUpdator() override = default;

  template <typename T>
  static void handleDispersiveMaterialADEUpdator(
      std::unordered_map<std::size_t, std::size_t>& map,
      std::vector<std::shared_ptr<T>>& dispersion_arr,
      const std::vector<std::shared_ptr<Material>>& material_arr) {
    std::size_t m_index = 0;
    std::size_t temp_i = 0;
    for (const auto& m : material_arr) {
      if (!m->dispersion()) {
        ++m_index;
        continue;
      }

      auto dispersive_material = std::dynamic_pointer_cast<T>(m);
      if (dispersive_material == nullptr) {
        ++m_index;
        continue;
      }

      map.insert({m_index, temp_i});
      dispersion_arr.emplace_back(dispersive_material);
      ++m_index;
      ++temp_i;
    }
  }
};

class LorentzADEUpdator : public LinearDispersiveMaterialADEUpdator {
 public:
  LorentzADEUpdator(std::shared_ptr<const GridSpace> grid_space,
                    std::shared_ptr<const CalculationParam> calculation_param,
                    std::shared_ptr<EMF> emf, IndexTask task);

  LorentzADEUpdator(const LorentzADEUpdator&) = delete;

  LorentzADEUpdator(LorentzADEUpdator&&) noexcept = default;

  LorentzADEUpdator& operator=(const LorentzADEUpdator&) = delete;

  LorentzADEUpdator& operator=(LorentzADEUpdator&&) noexcept = delete;

  ~LorentzADEUpdator() override = default;

  void updateE() override;

 protected:
 private:
  struct LorentzADERecord {
    Array4D<Real> _jx_prev, _jy_prev, _jz_prev, _jx, _jy, _jz;
  };
  std::unordered_map<std::size_t, std::size_t> _lorentz_map;
  std::vector<std::shared_ptr<LorentzMedium>> _lorentz_mediums;
  std::vector<ade::LorentzCoeff> _coeff;
  std::vector<LorentzADERecord> _j_record;
  Array3D<Real> _ex_prev, _ey_prev, _ez_prev;

  void init();

  auto updateEEdge() -> void override;
};

class DrudeADEUpdator : public LinearDispersiveMaterialADEUpdator {
 public:
  DrudeADEUpdator(std::shared_ptr<const GridSpace> grid_space,
                  std::shared_ptr<const CalculationParam> calculation_param,
                  std::shared_ptr<EMF> emf, IndexTask task);

  DrudeADEUpdator(const DrudeADEUpdator&) = delete;

  DrudeADEUpdator(DrudeADEUpdator&&) noexcept = default;

  DrudeADEUpdator& operator=(const DrudeADEUpdator&) = delete;

  DrudeADEUpdator& operator=(DrudeADEUpdator&&) noexcept = default;

  ~DrudeADEUpdator() override = default;

  void updateE() override;

 private:
  struct DrudeADERecord {
    Array4D<Real> _jx, _jy, _jz;
  };

  std::unordered_map<std::size_t, std::size_t> _drude_map;
  std::vector<std::shared_ptr<DrudeMedium>> _drude_mediums;
  std::vector<ade::DrudeCoeff> _coeff;
  std::vector<DrudeADERecord> _j_record;

  void init();

  auto updateEEdge() -> void override;
};

class DebyeADEUpdator : public LinearDispersiveMaterialADEUpdator {
 public:
  DebyeADEUpdator(std::shared_ptr<const GridSpace> grid_space,
                  std::shared_ptr<const CalculationParam> calculation_param,
                  std::shared_ptr<EMF> emf, IndexTask task);

  DebyeADEUpdator(const DebyeADEUpdator&) = delete;

  DebyeADEUpdator(DebyeADEUpdator&&) noexcept = default;

  DebyeADEUpdator& operator=(const DebyeADEUpdator&) = delete;

  DebyeADEUpdator& operator=(DebyeADEUpdator&&) noexcept = default;

  ~DebyeADEUpdator() override = default;

  void updateE() override;

 private:
  struct DebyeADERecord {
    Array4D<Real> _jx, _jy, _jz;
  };

  std::unordered_map<std::size_t, std::size_t> _debye_map;
  std::vector<std::shared_ptr<DebyeMedium>> _debye_mediums;
  std::vector<ade::DebyCoeff> _coeff;
  std::vector<DebyeADERecord> _j_record;

  void init();

  auto updateEEdge() -> void override;
};

}  // namespace xfdtd

#endif  // _XFDTD_CORE_DISPERSIVE_MATERIAL_UPDATOR_H_
