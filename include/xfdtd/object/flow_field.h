#ifndef __XFDTD_CORE_FLOW_FIELD_H__
#define __XFDTD_CORE_FLOW_FIELD_H__

#include <xfdtd/object/object.h>

#include <memory>
#include <string_view>
#include <vector>

#include "xfdtd/common/type_define.h"
#include "xfdtd/coordinate_system/coordinate_system.h"

namespace xfdtd {

class ModelShape;

class XFDTDFlowFieldException : public XFDTDException {
 public:
  explicit XFDTDFlowFieldException(
      std::string message = "XFDTD Flow Field Exception")
      : XFDTDException{std::move(message)} {}
};

class FlowFieldEntry {
 public:
  FlowFieldEntry() = default;
  FlowFieldEntry(Vector position, Real omega_p, Real gamma)
      : _position{position}, _omega_p{omega_p}, _gamma{gamma} {}

  auto position() const { return _position; }

  auto omegaP() const { return _omega_p; }

  auto gamma() const { return _gamma; }

 private:
  Vector _position{};
  Real _omega_p{}, _gamma{};
};

/**
 * @brief Read flow field
 * The format of the flow field is as follows:
 *  TITLE = "TEC3DS from TOPL3D at NT=10000, TAU=0.00000"
 *  VARIABLES = "X", "Y", "Z", "Fp", "Muc"
 *  ZONE T="Zone_1", I=70, J=108, K=60,
 *  F=POINT
 *    2.003403  0.2990603  0  4.773956e+10  8.218315e+09
 *    ...
 *
 */
class FlowField : public Object {
 public:
  FlowField(std::string_view name, std::string_view flow_field_model_info_path);

  FlowField(const FlowField&) = delete;

  FlowField(FlowField&&) noexcept = delete;

  FlowField& operator=(const FlowField&) = delete;

  FlowField& operator=(FlowField&&) noexcept = delete;

  virtual ~FlowField();

  std::string toString() const override;

  void init(std::shared_ptr<const GridSpace> grid_space,
            std::shared_ptr<CalculationParam> calculation_param,
            std::shared_ptr<EMF> emf) override;

  void correctMaterialSpace(std::size_t index) override;

  void correctUpdateCoefficient() override;

  auto flowFieldModelInfoPath() const -> std::string;

  auto objectModelInfoPath() const -> std::string;

  auto setFlowFieldModelInfoPath(std::string_view path) -> void;

  auto setObjectModelInfoPath(std::string_view path) -> void;

 private:
  std::string _flow_field_model_info_path{};
  std::vector<FlowFieldEntry> _flow_field_entries{};

  auto readModelInfo(std::string_view path) -> std::unique_ptr<ModelShape>;

  auto isInsideFlowField(const Vector& vector) const -> bool;

  auto isInsideObject(const Vector& vector) const -> bool;

  auto correctCoefficient(Index i, Index j, Index k,
                          const FlowFieldEntry& entry) -> void;

  struct DrudeMaterialParameter {
    Real epsilon_inf, omega_p, gamma;
  };

  auto getDrudeMaterialParameter(Index i, Index j,
                                 Index k) -> DrudeMaterialParameter;
};

}  // namespace xfdtd

#endif  // __XFDTD_CORE_FLOW_FIELD_H__
