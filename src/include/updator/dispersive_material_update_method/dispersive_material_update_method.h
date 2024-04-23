#ifndef __XFDTD_CORE_DISPERSIVE_MATERIAL_UPDATE_METHOD_H__
#define __XFDTD_CORE_DISPERSIVE_MATERIAL_UPDATE_METHOD_H__

#include <xfdtd/calculation_param/calculation_param.h>
#include <xfdtd/common/index_task.h>
#include <xfdtd/common/type_define.h>
#include <xfdtd/grid_space/grid_space.h>

namespace xfdtd {

class LinearDispersiveMaterialUpdateMethod {
 public:
  explicit LinearDispersiveMaterialUpdateMethod(Real epsilon_inf)
      : _epsilon_inf{epsilon_inf} {}

  virtual ~LinearDispersiveMaterialUpdateMethod() = default;

  virtual auto clone() const
      -> std::unique_ptr<LinearDispersiveMaterialUpdateMethod> = 0;

  virtual auto init(Real dt) -> void = 0;

  virtual auto correctCoeff(Index i, Index j, Index k,
                            const GridSpace* grid_space,
                            CalculationParam* calculation_param) const
      -> void = 0;

  virtual auto initUpdate(
      const GridSpace* grid_space,
      std::shared_ptr<const CalculationParam> calculation_param,
      std::shared_ptr<EMF> emf, Index m_index, const IndexTask& task)
      -> void = 0;

  virtual auto updateEx(Index i, Index j, Index k) -> void = 0;

  virtual auto updateEy(Index i, Index j, Index k) -> void = 0;

  virtual auto updateEz(Index i, Index j, Index k) -> void = 0;

  virtual auto updateTEM(Index i, Index j, Index k) -> void = 0;

 protected:
  Real _epsilon_inf;
  std::shared_ptr<const CalculationParam> _calculation_param;
  std::shared_ptr<EMF> _emf;
  Index _is, _js, _ks;

 private:
};

}  // namespace xfdtd

#endif  // __XFDTD_CORE_DISPERSIVE_MATERIAL_UPDATE_METHOD_H__
