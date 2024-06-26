#include <xfdtd/common/constant.h>
#include <xfdtd/common/type_define.h>
#include <xfdtd/coordinate_system/coordinate_system.h>
#include <xfdtd/electromagnetic_field/electromagnetic_field.h>
#include <xfdtd/grid_space/grid_space.h>
#include <xfdtd/nffft/nffft_time_domain.h>
#include <xfdtd/parallel/mpi_support.h>

#include <algorithm>
#include <cmath>
#include <filesystem>
#include <xtensor.hpp>
#include <xtensor/xnpy.hpp>

#include "nffft/nffft_td_data.h"

namespace xfdtd {

NFFFTTimeDomain::NFFFTTimeDomain(Index distance_x, Index distance_y,
                                 Index distance_z, Real theta, Real phi)
    : NFFFT{distance_x, distance_y, distance_z}, _theta{theta}, _phi{phi} {}

NFFFTTimeDomain::NFFFTTimeDomain(Index distance_x, Index distance_y,
                                 Index distance_z, std::string_view output_dir,
                                 Real theta, Real phi)
    : NFFFT{distance_x, distance_y, distance_z, output_dir},
      _theta{theta},
      _phi{phi} {}

NFFFTTimeDomain::~NFFFTTimeDomain() = default;

auto NFFFTTimeDomain::init(
    std::shared_ptr<const GridSpace> grid_space,
    std::shared_ptr<const CalculationParam> calculation_param,
    std::shared_ptr<const EMF> emf) -> void {
  defaultInit(grid_space, calculation_param, emf);
  _td_plane_xn = std::make_shared<TDPlaneData<Axis::Direction::XN>>(
      grid_space, calculation_param, emf, nodeTaskSurfaceXN());
  _td_plane_xp = std::make_shared<TDPlaneData<Axis::Direction::XP>>(
      grid_space, calculation_param, emf, nodeTaskSurfaceXP());
  _td_plane_yn = std::make_shared<TDPlaneData<Axis::Direction::YN>>(
      grid_space, calculation_param, emf, nodeTaskSurfaceYN());
  _td_plane_yp = std::make_shared<TDPlaneData<Axis::Direction::YP>>(
      grid_space, calculation_param, emf, nodeTaskSurfaceYP());
  _td_plane_zn = std::make_shared<TDPlaneData<Axis::Direction::ZN>>(
      grid_space, calculation_param, emf, nodeTaskSurfaceZN());
  _td_plane_zp = std::make_shared<TDPlaneData<Axis::Direction::ZP>>(
      grid_space, calculation_param, emf, nodeTaskSurfaceZP());
}

auto NFFFTTimeDomain::initTimeDependentVariable() -> void {
  auto distance_range_e = distanceRange<EMF::Attribute::E>();
  auto distance_range_h = distanceRange<EMF::Attribute::H>();

  auto aux_e = (distance_range_e.size()) /
               (constant::C_0 * calculationParam()->timeParam()->dt());
  auto aux_h = (distance_range_h.size()) /
               (constant::C_0 * calculationParam()->timeParam()->dt());

  auto num_e = static_cast<Index>(std::ceil(aux_e)) +
               calculationParam()->timeParam()->size() + 4;  // why +4?
  auto num_h = static_cast<Index>(std::ceil(aux_h)) +
               calculationParam()->timeParam()->size() + 4;

  auto observation_direction = observationDirection();

  _td_plane_xn->initForUpdate(num_e, num_h, distance_range_e, distance_range_h,
                              observation_direction);
  _td_plane_xp->initForUpdate(num_e, num_h, distance_range_e, distance_range_h,
                              observation_direction);
  _td_plane_yn->initForUpdate(num_e, num_h, distance_range_e, distance_range_h,
                              observation_direction);
  _td_plane_yp->initForUpdate(num_e, num_h, distance_range_e, distance_range_h,
                              observation_direction);
  _td_plane_zn->initForUpdate(num_e, num_h, distance_range_e, distance_range_h,
                              observation_direction);
  _td_plane_zp->initForUpdate(num_e, num_h, distance_range_e, distance_range_h,
                              observation_direction);
}

auto NFFFTTimeDomain::update() -> void {
  _td_plane_xn->update<EMF::Attribute::E>();
  _td_plane_xp->update<EMF::Attribute::E>();
  _td_plane_yn->update<EMF::Attribute::E>();
  _td_plane_yp->update<EMF::Attribute::E>();
  _td_plane_zn->update<EMF::Attribute::E>();
  _td_plane_zp->update<EMF::Attribute::E>();

  _td_plane_xn->update<EMF::Attribute::H>();
  _td_plane_xp->update<EMF::Attribute::H>();
  _td_plane_yn->update<EMF::Attribute::H>();
  _td_plane_yp->update<EMF::Attribute::H>();
  _td_plane_zn->update<EMF::Attribute::H>();
  _td_plane_zp->update<EMF::Attribute::H>();
}

auto NFFFTTimeDomain::processFarField() const -> void {
  if (!valid()) {
    return;
  }

  auto node_w_theta = wx() * cosTheta() * cosPhi() +
                      wy() * cosTheta() * sinPhi() - wz() * sinTheta();
  auto node_w_phi = -wx() * sinPhi() + wy() * cosPhi();

  auto node_u_theta = ux() * cosTheta() * cosPhi() +
                      uy() * cosTheta() * sinPhi() - uz() * sinTheta();
  auto node_u_phi = -ux() * sinPhi() + uy() * cosPhi();

  Array1D<Real> node_data;
  node_data = xt::concatenate(xt::xtuple(node_data, node_w_theta));
  node_data = xt::concatenate(xt::xtuple(node_data, node_w_phi));
  node_data = xt::concatenate(xt::xtuple(node_data, node_u_theta));
  node_data = xt::concatenate(xt::xtuple(node_data, node_u_phi));

  auto nffft_gather_func = [this](const auto& send_data, auto&& recv_data) {
    if (this->nffftMPIConfig().size() <= 1) {
      recv_data = send_data;
      return;
    }

    MpiSupport::instance().reduceSum(this->nffftMPIConfig(), send_data.data(),
                                     recv_data.data(), send_data.size());
  };

  auto data = xt::zeros_like(node_data);
  nffft_gather_func(node_data, data);

  if (!nffftMPIConfig().isRoot()) {
    return;
  }

  const auto output_dir = outputDir();

  if (!std::filesystem::exists(output_dir)) {
    std::filesystem::create_directories(output_dir);
  }

  auto&& w_theta = xt::view(data, xt::range(0, node_w_theta.size()));
  auto&& w_phi = xt::view(
      data,
      xt::range(node_w_theta.size(), node_w_theta.size() + node_w_phi.size()));
  auto&& u_theta = xt::view(
      data,
      xt::range(node_w_theta.size() + node_w_phi.size(),
                node_w_theta.size() + node_w_phi.size() + node_u_theta.size()));
  auto&& u_phi = xt::view(
      data,
      xt::range(node_w_theta.size() + node_w_phi.size() + node_u_theta.size(),
                _));

  xt::dump_npy((std::filesystem::path{output_dir} / "w_theta.npy").string(),
               w_theta);
  xt::dump_npy((std::filesystem::path{output_dir} / "w_phi.npy").string(),
               w_phi);
  xt::dump_npy((std::filesystem::path{output_dir} / "u_theta.npy").string(),
               u_theta);
  xt::dump_npy((std::filesystem::path{output_dir} / "u_phi.npy").string(),
               u_phi);
}

template <EMF::Attribute attribute>
auto NFFFTTimeDomain::distanceRange() const -> Range<Real> {
  auto distance_range = [this](const auto& observer_direction,
                               const auto& global_task, Axis::XYZ xyz) {
    auto is = global_task.xRange().start();
    auto js = global_task.yRange().start();
    auto ks = global_task.zRange().start();
    auto ie = global_task.xRange().end();
    auto je = global_task.yRange().end();
    auto ke = global_task.zRange().end();

    auto min_dis = std::numeric_limits<Real>::max();
    auto max_dis = std::numeric_limits<Real>::lowest();

    auto r_vector = [](Index i, Index j, Index k, const GridSpace* grid_space,
                       Axis::XYZ xyz) {
      if (xyz == Axis::XYZ::X) {
        return Vector{grid_space->eNodeX()(i), grid_space->hNodeY()(j),
                      grid_space->hNodeZ()(k)};
      }
      if (xyz == Axis::XYZ::Y) {
        return Vector{grid_space->hNodeX()(i), grid_space->eNodeY()(j),
                      grid_space->hNodeZ()(k)};
      }
      if (xyz == Axis::XYZ::Z) {
        return Vector{grid_space->hNodeX()(i), grid_space->hNodeY()(j),
                      grid_space->eNodeZ()(k)};
      }
      throw std::runtime_error("Invalid Axis::XYZ");
    };

    for (auto k = ks; k < ke; ++k) {
      for (auto j = js; j < je; ++j) {
        for (auto i = is; i < ie; ++i) {
          auto r = r_vector(i, j, k, this->gridSpace()->globalGridSpace().get(),
                            xyz);
          auto dis = r.x() * observer_direction.x() +
                     r.y() * observer_direction.y() +
                     r.z() * observer_direction.z();
          min_dis = std::min(min_dis, dis);
          max_dis = std::max(max_dis, dis);
        }
      }
    }

    return Range<Real>{min_dis, max_dis};
  };

  auto observation_direction = observationDirection();

  auto range_xn = distance_range(observation_direction, globalTaskSurfaceXN(),
                                 Axis::XYZ::X);
  auto range_xp = distance_range(observation_direction, globalTaskSurfaceXP(),
                                 Axis::XYZ::X);
  auto range_yn = distance_range(observation_direction, globalTaskSurfaceYN(),
                                 Axis::XYZ::Y);
  auto range_yp = distance_range(observation_direction, globalTaskSurfaceYP(),
                                 Axis::XYZ::Y);
  auto range_zn = distance_range(observation_direction, globalTaskSurfaceZN(),
                                 Axis::XYZ::Z);
  auto range_zp = distance_range(observation_direction, globalTaskSurfaceZP(),
                                 Axis::XYZ::Z);

  return {std::min({range_xn.start(), range_xp.start(), range_yn.start(),
                    range_yp.start(), range_zn.start(), range_zp.start()}),
          std::max({range_xn.end(), range_xp.end(), range_yn.end(),
                    range_yp.end(), range_zn.end(), range_zp.end()})};
}

auto NFFFTTimeDomain::wx() const -> Array1D<Real> {
  return _td_plane_zn->wa() + _td_plane_zp->wa() + _td_plane_yn->wb() +
         _td_plane_yp->wb();
}

auto NFFFTTimeDomain::wy() const -> Array1D<Real> {
  return _td_plane_xn->wa() + _td_plane_xp->wa() + _td_plane_zn->wb() +
         _td_plane_zp->wb();
}

auto NFFFTTimeDomain::wz() const -> Array1D<Real> {
  return _td_plane_xn->wb() + _td_plane_xp->wb() + _td_plane_yn->wa() +
         _td_plane_yp->wa();
}

auto NFFFTTimeDomain::ux() const -> Array1D<Real> {
  return _td_plane_zn->ua() + _td_plane_zp->ua() + _td_plane_yn->ub() +
         _td_plane_yp->ub();
}

auto NFFFTTimeDomain::uy() const -> Array1D<Real> {
  return _td_plane_xn->ua() + _td_plane_xp->ua() + _td_plane_zn->ub() +
         _td_plane_zp->ub();
}

auto NFFFTTimeDomain::uz() const -> Array1D<Real> {
  return _td_plane_xn->ub() + _td_plane_xp->ub() + _td_plane_yn->ua() +
         _td_plane_yp->ua();
}

auto NFFFTTimeDomain::observationDirection() const -> Vector {
  return Vector{sinTheta() * cosPhi(), sinTheta() * sinPhi(), cosTheta()};
}

auto NFFFTTimeDomain::cosTheta() const -> Real { return std::cos(_theta); }

auto NFFFTTimeDomain::sinTheta() const -> Real { return std::sin(_theta); }

auto NFFFTTimeDomain::cosPhi() const -> Real { return std::cos(_phi); }

auto NFFFTTimeDomain::sinPhi() const -> Real { return std::sin(_phi); }

}  // namespace xfdtd