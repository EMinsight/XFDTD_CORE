// #ifndef _XFDTD_CORE_DISPERSIVE_MATERIAL_UPDATOR_H_
// #define _XFDTD_CORE_DISPERSIVE_MATERIAL_UPDATOR_H_

// #include <updator/basic_updator.h>
// #include <xfdtd/calculation_param/calculation_param.h>
// #include <xfdtd/common/index_task.h>
// #include <xfdtd/common/type_define.h>
// #include <xfdtd/material/dispersive_material.h>

// #include <memory>
// #include <vector>

// namespace xfdtd {

// class XFDTDLinerDispersiveMaterialUpdatorException : public XFDTDException {
//  public:
//   explicit XFDTDLinerDispersiveMaterialUpdatorException(
//       const std::string& message)
//       : XFDTDException(message) {}
// };

// class LinearDispersiveMaterialUpdateMethod;

// class LinearDispersiveMaterialUpdator : public BasicUpdator3D {
//  public:
//  public:
//   LinearDispersiveMaterialUpdator(
//       std::vector<std::shared_ptr<Material>> material_arr,
//       std::shared_ptr<const GridSpace> grid_space,
//       std::shared_ptr<const CalculationParam> calculation_param,
//       std::shared_ptr<EMF> emf, IndexTask task);

//   ~LinearDispersiveMaterialUpdator() override = default;

//   auto updateE() -> void override;

//  private:
//   std::vector<Index> _map;
//   std::vector<std::shared_ptr<LinearDispersiveMaterialUpdateMethod>>
//       _update_methods;

//   auto init(std::vector<std::shared_ptr<Material>> material_arr) -> void;
// };

// class LinearDispersiveMaterialUpdator1D : public BasicUpdatorTEM {
//  public:
//   LinearDispersiveMaterialUpdator1D(
//       const std::vector<std::shared_ptr<Material>>& material_arr,
//       std::shared_ptr<const GridSpace> grid_space,
//       std::shared_ptr<const CalculationParam> calculation_param,
//       std::shared_ptr<EMF> emf, IndexTask task);

//   ~LinearDispersiveMaterialUpdator1D() override = default;

//   auto updateE() -> void override;

//  private:
//   std::vector<Index> _map;
//   std::vector<std::shared_ptr<LinearDispersiveMaterialUpdateMethod>>
//       _update_methods;

//   auto init(const std::vector<std::shared_ptr<Material>>& material_arr) -> void;
// };

// }  // namespace xfdtd

// #endif  // _XFDTD_CORE_DISPERSIVE_MATERIAL_UPDATOR_H_
