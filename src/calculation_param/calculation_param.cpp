#include <xfdtd/calculation_param/calculation_param.h>
#include <xfdtd/common/constant.h>
#include <xfdtd/grid_space/grid_space.h>

#include <xtensor.hpp>

namespace xfdtd {

XFDTDCalculationParamException::XFDTDCalculationParamException(
    std::string message)
    : XFDTDException{std::move(message)} {}

CalculationParam::CalculationParam()
    : _time_param{nullptr},
      _material_param{nullptr},
      _fdtd_coefficient{std::make_unique<FDTDUpdateCoefficient>()} {}

const std::unique_ptr<TimeParam>& CalculationParam::timeParam() const {
  return _time_param;
}

const std::unique_ptr<MaterialParam>& CalculationParam::materialParam() const {
  return _material_param;
}

const std::unique_ptr<FDTDUpdateCoefficient>&
CalculationParam::fdtdCoefficient() const {
  return _fdtd_coefficient;
}

std::unique_ptr<TimeParam>& CalculationParam::timeParam() {
  return _time_param;
}

std::unique_ptr<MaterialParam>& CalculationParam::materialParam() {
  return _material_param;
}

std::unique_ptr<FDTDUpdateCoefficient>& CalculationParam::fdtdCoefficient() {
  return _fdtd_coefficient;
}

void CalculationParam::generateMaterialSpaceParam(const GridSpace* grid_space) {
  // TODO(franzero): don't support to nonuniform grid space
  if (grid_space->type() != GridSpace::Type::UNIFORM) {
    throw XFDTDCalculationParamException{
        "Nonuniform grid space is not "
        "supported yet"};
  }

  auto nx{grid_space->sizeX()};
  auto ny{grid_space->sizeY()};
  auto nz{grid_space->sizeZ()};
  auto& eps_x{materialParam()->epsX()};
  xt::view(eps_x, xt::range(0, nx), xt::range(1, ny), xt::range(1, nz)) =
      0.25 *
      (xt::view(eps_x, xt::range(0, nx), xt::range(1, ny), xt::range(1, nz)) +
       xt::view(eps_x, xt::range(0, nx), xt::range(0, ny - 1),
                xt::range(1, nz)) +
       xt::view(eps_x, xt::range(0, nx), xt::range(1, ny),
                xt::range(0, nz - 1)) +
       xt::view(eps_x, xt::range(0, nx), xt::range(0, ny - 1),
                xt::range(0, nz - 1)));
  auto& eps_y{materialParam()->epsY()};
  xt::view(eps_y, xt::range(1, nx), xt::range(0, ny), xt::range(1, nz)) =
      0.25 *
      (xt::view(eps_y, xt::range(1, nx), xt::range(0, ny), xt::range(1, nz)) +
       xt::view(eps_y, xt::range(0, nx - 1), xt::range(0, ny),
                xt::range(1, nz)) +
       xt::view(eps_y, xt::range(1, nx), xt::range(0, ny),
                xt::range(0, nz - 1)) +
       xt::view(eps_y, xt::range(0, nx - 1), xt::range(0, ny),
                xt::range(0, nz - 1)));
  auto& eps_z{materialParam()->epsZ()};
  xt::view(eps_z, xt::range(1, nx), xt::range(1, ny), xt::range(0, nz)) =
      0.25 *
      (xt::view(eps_z, xt::range(1, nx), xt::range(1, ny), xt::range(0, nz)) +
       xt::view(eps_z, xt::range(0, nx - 1), xt::range(1, ny),
                xt::range(0, nz)) +
       xt::view(eps_z, xt::range(1, nx), xt::range(0, ny - 1),
                xt::range(0, nz)) +
       xt::view(eps_z, xt::range(0, nx - 1), xt::range(0, ny - 1),
                xt::range(0, nz)));
  auto& mu_x{materialParam()->muX()};
  xt::view(mu_x, xt::range(1, nx), xt::range(0, ny), xt::range(0, nz)) =
      2 * xt::view(mu_x, xt::range(1, nx), xt::range(0, ny), xt::range(0, nz)) *
      xt::view(mu_x, xt::range(0, nx - 1), xt::range(0, ny), xt::range(0, nz)) /
      (xt::view(mu_x, xt::range(1, nx), xt::range(0, ny), xt::range(0, nz)) +
       xt::view(mu_x, xt::range(0, nx - 1), xt::range(0, ny),
                xt::range(0, nz)));
  auto& mu_y{materialParam()->muY()};
  xt::view(mu_y, xt::range(0, nx), xt::range(1, ny), xt::range(0, nz)) =
      2 * xt::view(mu_y, xt::range(0, nx), xt::range(1, ny), xt::range(0, nz)) *
      xt::view(mu_y, xt::range(0, nx), xt::range(0, ny - 1), xt::range(0, nz)) /
      (xt::view(mu_y, xt::range(0, nx), xt::range(1, ny), xt::range(0, nz)) +
       xt::view(mu_y, xt::range(0, nx), xt::range(0, ny - 1),
                xt::range(0, nz)));
  auto& mu_z{materialParam()->muZ()};
  xt::view(mu_z, xt::range(0, nx), xt::range(0, ny), xt::range(1, nz)) =
      2 * xt::view(mu_z, xt::range(0, nx), xt::range(0, ny), xt::range(1, nz)) *
      xt::view(mu_z, xt::range(0, nx), xt::range(0, ny), xt::range(0, nz - 1)) /
      (xt::view(mu_z, xt::range(0, nx), xt::range(0, ny), xt::range(1, nz)) +
       xt::view(mu_z, xt::range(0, nx), xt::range(0, ny),
                xt::range(0, nz - 1)));
  auto& sigma_e_x{materialParam()->sigmaEX()};
  xt::view(sigma_e_x, xt::range(0, nx), xt::range(1, ny), xt::range(1, nz)) =
      0.25 * (xt::view(sigma_e_x, xt::range(0, nx), xt::range(1, ny),
                       xt::range(1, nz)) +
              xt::view(sigma_e_x, xt::range(0, nx), xt::range(0, ny - 1),
                       xt::range(1, nz)) +
              xt::view(sigma_e_x, xt::range(0, nx), xt::range(1, ny),
                       xt::range(0, nz - 1)) +
              xt::view(sigma_e_x, xt::range(0, nx), xt::range(0, ny - 1),
                       xt::range(0, nz - 1)));
  auto& sigma_e_y{materialParam()->sigmaEY()};
  xt::view(sigma_e_y, xt::range(1, nx), xt::range(0, ny), xt::range(1, nz)) =
      0.25 * (xt::view(sigma_e_y, xt::range(1, nx), xt::range(0, ny),
                       xt::range(1, nz)) +
              xt::view(sigma_e_y, xt::range(0, nx - 1), xt::range(0, ny),
                       xt::range(1, nz)) +
              xt::view(sigma_e_y, xt::range(1, nx), xt::range(0, ny),
                       xt::range(0, nz - 1)) +
              xt::view(sigma_e_y, xt::range(0, nx - 1), xt::range(0, ny),
                       xt::range(0, nz - 1)));
  auto& sigma_e_z{materialParam()->sigmaEZ()};
  xt::view(sigma_e_z, xt::range(1, nx), xt::range(1, ny), xt::range(0, nz)) =
      0.25 * (xt::view(sigma_e_z, xt::range(1, nx), xt::range(1, ny),
                       xt::range(0, nz)) +
              xt::view(sigma_e_z, xt::range(0, nx - 1), xt::range(1, ny),
                       xt::range(0, nz)) +
              xt::view(sigma_e_z, xt::range(1, nx), xt::range(0, ny - 1),
                       xt::range(0, nz)) +
              xt::view(sigma_e_z, xt::range(0, nx - 1), xt::range(0, ny - 1),
                       xt::range(0, nz)));
  auto& sigma_m_x{materialParam()->sigmaMX()};
  xt::view(sigma_m_x, xt::range(1, nx), xt::range(0, ny), xt::range(0, nz)) =
      2 *
      xt::view(sigma_m_x, xt::range(1, nx), xt::range(0, ny),
               xt::range(0, nz)) *
      xt::view(sigma_m_x, xt::range(0, nx - 1), xt::range(0, ny),
               xt::range(0, nz)) /
      (xt::view(sigma_m_x, xt::range(1, nx), xt::range(0, ny),
                xt::range(0, nz)) +
       xt::view(sigma_m_x, xt::range(0, nx - 1), xt::range(0, ny),
                xt::range(0, nz)));
  auto& sigma_m_y{materialParam()->sigmaMY()};
  xt::view(sigma_m_y, xt::range(0, nx), xt::range(1, ny), xt::range(0, nz)) =
      2 *
      xt::view(sigma_m_y, xt::range(0, nx), xt::range(1, ny),
               xt::range(0, nz)) *
      xt::view(sigma_m_y, xt::range(0, nx), xt::range(0, ny - 1),
               xt::range(0, nz)) /
      (xt::view(sigma_m_y, xt::range(0, nx), xt::range(1, ny),
                xt::range(0, nz)) +
       xt::view(sigma_m_y, xt::range(0, nx), xt::range(0, ny - 1),
                xt::range(0, nz)));
  auto& sigma_m_z{materialParam()->sigmaMZ()};
  xt::view(sigma_m_z, xt::range(0, nx), xt::range(0, ny), xt::range(1, nz)) =
      2 *
      xt::view(sigma_m_z, xt::range(0, nx), xt::range(0, ny),
               xt::range(1, nz)) *
      xt::view(sigma_m_z, xt::range(0, nx), xt::range(0, ny),
               xt::range(0, nz - 1)) /
      (xt::view(sigma_m_z, xt::range(0, nx), xt::range(0, ny),
                xt::range(1, nz)) +
       xt::view(sigma_m_z, xt::range(0, nx), xt::range(0, ny),
                xt::range(0, nz - 1)));
}

void CalculationParam::calculateCoefficient(const GridSpace* grid_space) {
  const auto& h_size_x{grid_space->hSizeX()};
  const auto& h_size_y{grid_space->hSizeY()};
  const auto& h_size_z{grid_space->hSizeZ()};
  const auto& e_size_x{grid_space->eSizeX()};
  const auto& e_size_y{grid_space->eSizeY()};
  const auto& e_size_z{grid_space->eSizeZ()};

  auto get_xyz = [](const auto& size_x, const auto& size_y,
                    const auto& size_z) {
    return xt::meshgrid(size_x, size_y, size_z);
  };

  auto calculate_coff_e{[](const auto& da, const auto& db, const auto& dc,
                           const auto& dt, const auto& eps, const auto& sigma,
                           auto&& cece, auto&& cecha, auto&& cechb) {
    cece = (2 * eps - dt * sigma) / (2 * eps + dt * sigma);
    cecha = -(2 * dt / db) / (2 * eps + sigma * dt);
    cechb = (2 * dt / da) / (2 * eps + sigma * dt);
  }};

  auto calculate_coff_h{[](const auto& da, const auto& db, const auto& dc,
                           const auto& dt, const auto& mu, const auto& sigma_m,
                           auto&& chch, auto&& chcea, auto&& chceb) {
    chch = (2 * mu - dt * sigma_m) / (2 * mu + dt * sigma_m);
    chcea = (2 * dt / db) / (2 * mu + sigma_m * dt);
    chceb = -(2 * dt / da) / (2 * mu + sigma_m * dt);
  }};

  auto dt{timeParam()->dt()};
  auto nx{grid_space->sizeX()};
  auto ny{grid_space->sizeY()};
  auto nz{grid_space->sizeZ()};
  auto& eps_x{materialParam()->epsX()};
  auto& eps_y{materialParam()->epsY()};
  auto& eps_z{materialParam()->epsZ()};
  auto& mu_x{materialParam()->muX()};
  auto& mu_y{materialParam()->muY()};
  auto& mu_z{materialParam()->muZ()};
  auto& sigma_e_x{materialParam()->sigmaEX()};
  auto& sigma_e_y{materialParam()->sigmaEY()};
  auto& sigma_e_z{materialParam()->sigmaEZ()};
  auto& sigma_m_x{materialParam()->sigmaMX()};
  auto& sigma_m_y{materialParam()->sigmaMY()};
  auto& sigma_m_z{materialParam()->sigmaMZ()};

  auto& cexe{fdtdCoefficient()->cexe()};
  auto& cexhy{fdtdCoefficient()->cexhy()};
  auto& cexhz{fdtdCoefficient()->cexhz()};
  auto& ceye{fdtdCoefficient()->ceye()};
  auto& ceyhx{fdtdCoefficient()->ceyhx()};
  auto& ceyhz{fdtdCoefficient()->ceyhz()};
  auto& ceze{fdtdCoefficient()->ceze()};
  auto& cezhx{fdtdCoefficient()->cezhx()};
  auto& cezhy{fdtdCoefficient()->cezhy()};
  auto& chxh{fdtdCoefficient()->chxh()};
  auto& chxey{fdtdCoefficient()->chxey()};
  auto& chxez{fdtdCoefficient()->chxez()};
  auto& chyh{fdtdCoefficient()->chyh()};
  auto& chyez{fdtdCoefficient()->chyez()};
  auto& chyex{fdtdCoefficient()->chyex()};
  auto& chzh{fdtdCoefficient()->chzh()};
  auto& chzex{fdtdCoefficient()->chzex()};
  auto& chzey{fdtdCoefficient()->chzey()};

  auto e_x_size{get_xyz(e_size_x, h_size_y, h_size_z)};
  calculate_coff_e(std::get<1>(e_x_size), std::get<2>(e_x_size),
                   std::get<0>(e_x_size), dt, eps_x, sigma_e_x, cexe, cexhy,
                   cexhz);

  auto e_y_size{get_xyz(h_size_x, e_size_y, h_size_z)};
  calculate_coff_e(std::get<2>(e_y_size), std::get<0>(e_y_size),
                   std::get<1>(e_y_size), dt, eps_y, sigma_e_y, ceye, ceyhz,
                   ceyhx);

  auto e_z_size{get_xyz(h_size_x, h_size_y, e_size_z)};
  calculate_coff_e(std::get<0>(e_z_size), std::get<1>(e_z_size),
                   std::get<2>(e_z_size), dt, eps_z, sigma_e_z, ceze, cezhx,
                   cezhy);

  auto h_x_size{get_xyz(h_size_x, e_size_y, e_size_z)};
  calculate_coff_h(std::get<1>(h_x_size), std::get<2>(h_x_size),
                   std::get<0>(h_x_size), dt, mu_x, sigma_m_x, chxh, chxey,
                   chxez);

  auto h_y_size{get_xyz(e_size_x, h_size_y, e_size_z)};
  calculate_coff_h(std::get<2>(h_y_size), std::get<0>(h_y_size),
                   std::get<1>(h_y_size), dt, mu_y, sigma_m_y, chyh, chyez,
                   chyex);

  auto h_z_size{get_xyz(e_size_x, e_size_y, h_size_z)};
  calculate_coff_h(std::get<0>(h_z_size), std::get<1>(h_z_size),
                   std::get<2>(h_z_size), dt, mu_z, sigma_m_z, chzh, chzex,
                   chzey);
}

void CalculationParam::setTimeParam(std::unique_ptr<TimeParam> time_param) {
  _time_param = std::move(time_param);
}

void CalculationParam::setMaterialParam(
    std::unique_ptr<MaterialParam> material_param) {
  _material_param = std::move(material_param);
}

// void CalculationParam::allocate(const GridSpace* grid_space) {
//   auto nx{grid_space->hNodeX().size()};
//   auto ny{grid_space->hNodeY().size()};
//   auto nz{grid_space->hNodeZ().size()};

//   materialParam()->epsX().resize({nx, ny + 1, nz + 1});
//   materialParam()->epsY().resize({nx + 1, ny, nz + 1});
//   materialParam()->epsZ().resize({nx + 1, ny + 1, nz});
//   materialParam()->epsX().fill(constant::EPSILON_0);
//   materialParam()->epsY().fill(constant::EPSILON_0);
//   materialParam()->epsZ().fill(constant::EPSILON_0);

//   materialParam()->muX().resize({nx + 1, ny, nz});
//   materialParam()->muY().resize({nx, ny + 1, nz});
//   materialParam()->muZ().resize({nx, ny, nz + 1});
//   materialParam()->muX().fill(constant::MU_0);
//   materialParam()->muY().fill(constant::MU_0);
//   materialParam()->muZ().fill(constant::MU_0);

//   materialParam()->sigmaEX().resize({nx, ny + 1, nz + 1});
//   materialParam()->sigmaEY().resize({nx + 1, ny, nz + 1});
//   materialParam()->sigmaEZ().resize({nx + 1, ny + 1, nz});
//   materialParam()->sigmaEX().fill(1e-20);
//   materialParam()->sigmaEY().fill(1e-20);
//   materialParam()->sigmaEZ().fill(1e-20);

//   materialParam()->sigmaMX().resize({nx + 1, ny, nz});
//   materialParam()->sigmaMY().resize({nx, ny + 1, nz});
//   materialParam()->sigmaMZ().resize({nx, ny, nz + 1});
//   materialParam()->sigmaMX().fill(1e-20);
//   materialParam()->sigmaMY().fill(1e-20);
//   materialParam()->sigmaMZ().fill(1e-20);
// }

}  // namespace xfdtd
