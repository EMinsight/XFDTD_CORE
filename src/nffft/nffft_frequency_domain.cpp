#include <xfdtd/nffft/nffft_frequency_domain.h>
#include <xfdtd/parallel/mpi_support.h>

#include <filesystem>
#include <iomanip>
#include <xtensor.hpp>
#include <xtensor/xnpy.hpp>

#include "nffft/nffft_fd_data.h"

namespace xfdtd {

NFFFTFrequencyDomain::NFFFTFrequencyDomain(Index distance_x, Index distance_y,
                                           Index distance_z,
                                           Array1D<Real> frequencies)
    : NFFFT{distance_x, distance_y, distance_z},
      _frequencies{std::move(frequencies)} {}

NFFFTFrequencyDomain::NFFFTFrequencyDomain(Index distance_x, Index distance_y,
                                           Index distance_z,
                                           Array1D<Real> frequencies,
                                           std::string_view output_dir)
    : NFFFT{distance_x, distance_y, distance_z, output_dir},
      _frequencies{std::move(frequencies)} {}

NFFFTFrequencyDomain::~NFFFTFrequencyDomain() = default;

auto NFFFTFrequencyDomain::init(
    std::shared_ptr<const GridSpace> grid_space,
    std::shared_ptr<const CalculationParam> calculation_param,
    std::shared_ptr<const EMF> emf) -> void {
  defaultInit(grid_space, calculation_param, emf);
  generateSurface();
}

void NFFFTFrequencyDomain::initTimeDependentVariable() {
  const auto total_time_step{calculationParam()->timeParam()->endTimeStep() -
                             calculationParam()->timeParam()->startTimeStep()};
  const auto dt{calculationParam()->timeParam()->dt()};
  for (auto&& fd_data : _fd_plane_data) {
    fd_data.initDFT(total_time_step, dt);
  }
}

auto NFFFTFrequencyDomain::update() -> void {
  auto current_time_step = calculationParam()->timeParam()->currentTimeStep();

  for (auto f{0}; f < _frequencies.size(); ++f) {
    for (auto&& fd_data : _fd_plane_data) {
      fd_data.update(current_time_step);
    }
  }
}

auto NFFFTFrequencyDomain::processFarField(const Array1D<Real>& theta, Real phi,
                                           const std::string& sub_dir,
                                           const Vector& origin) const -> void {
  processFarField(theta, Array1D<Real>{phi}, sub_dir, origin);
}

auto NFFFTFrequencyDomain::processFarField(Real theta, const Array1D<Real>& phi,
                                           const std::string& sub_dir,
                                           const Vector& origin) const -> void {
  processFarField(Array1D<Real>{theta}, phi, sub_dir, origin);
}

void NFFFTFrequencyDomain::outputRadiationPower() {
  if (!valid()) {
    return;
  }

  auto num_freq = _fd_plane_data.size();
  Array1D<Real> freq_arr = xt::zeros<Real>({num_freq});
  Array1D<Real> node_power_arr = xt::zeros<Real>({num_freq});
  Array1D<Real> power_arr = xt::zeros<Real>({num_freq});

  for (auto i = 0; i < num_freq; ++i) {
    freq_arr(i) = _fd_plane_data[i].frequency();
    node_power_arr(i) = _fd_plane_data[i].power();
  }

  if (nffftMPIConfig().size() <= 1) {
    power_arr = node_power_arr;
  } else {
    MpiSupport::instance().reduceSum(nffftMPIConfig(), node_power_arr.data(),
                                     power_arr.data(), node_power_arr.size());
  }

  if (!nffftMPIConfig().isRoot()) {
    return;
  }

  const auto output_dir = outputDir();
  if (!std::filesystem::exists(output_dir)) {
    std::filesystem::create_directories(output_dir);
  }

  auto data = xt::stack(xt::xtuple(freq_arr, power_arr));
  xt::dump_npy((std::filesystem::path{output_dir} / "power.npy").string(),
               data);
}

auto NFFFTFrequencyDomain::processFarField(const Array1D<Real>& theta,
                                           const Array1D<Real>& phi,
                                           const std::string& sub_dir,
                                           const Vector& origin) const -> void {
  if (!valid()) {
    return;
  }

  Array1D<std::complex<Real>> node_data;
  for (const auto& f : _fd_plane_data) {
    const auto freq = f.frequency();
    auto node_a_theta = f.aTheta(theta, phi, origin);
    auto node_f_phi = f.fPhi(theta, phi, origin);
    auto node_a_phi = f.aPhi(theta, phi, origin);
    auto node_f_theta = f.fTheta(theta, phi, origin);

    Array1D<std::complex<Real>> a_theta = xt::zeros_like(node_a_theta);
    Array1D<std::complex<Real>> f_phi = xt::zeros_like(node_f_phi);
    Array1D<std::complex<Real>> a_phi = xt::zeros_like(node_a_phi);
    Array1D<std::complex<Real>> f_theta = xt::zeros_like(node_f_theta);

    node_data = xt::concatenate(xt::xtuple(node_data, node_a_theta));
    node_data = xt::concatenate(xt::xtuple(node_data, node_f_phi));
    node_data = xt::concatenate(xt::xtuple(node_data, node_a_phi));
    node_data = xt::concatenate(xt::xtuple(node_data, node_f_theta));
  }

  // rename later
  auto nffft_gather_func = [this](const Array1D<std::complex<Real>>& send_data,
                                  Array1D<std::complex<Real>>& recv_data) {
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

  auto output_dir = outputDir();
  const auto offsets = theta.size() * phi.size();
  for (auto i{0}; i < _fd_plane_data.size(); ++i) {
    auto freq = _fd_plane_data[i].frequency();
    std::stringstream ss;
    ss << std::fixed << std::setprecision(2);
    ss << freq / 1e9 << "GHz";
    const auto prefix_file_name = ss.str();
    if (!std::filesystem::exists(std::filesystem::path{output_dir} / sub_dir)) {
      std::filesystem::create_directories(std::filesystem::path{output_dir} /
                                          sub_dir);
    }

    auto&& fd_data =
        xt::view(data, xt::range(4 * i * offsets, (4 * i + 4) * offsets));
    auto a_theta = xt::view(fd_data, xt::range(0, offsets));
    auto f_phi = xt::view(fd_data, xt::range(offsets, 2 * offsets));
    auto a_phi = xt::view(fd_data, xt::range(2 * offsets, 3 * offsets));
    auto f_theta = xt::view(fd_data, xt::range(3 * offsets, 4 * offsets));

    xt::dump_npy((std::filesystem::path{output_dir} / sub_dir /
                  (prefix_file_name + "_a_theta.npy"))
                     .string(),
                 a_theta);
    xt::dump_npy((std::filesystem::path{output_dir} / sub_dir /
                  (prefix_file_name + "_f_phi.npy"))
                     .string(),
                 f_phi);
    xt::dump_npy((std::filesystem::path{output_dir} / sub_dir /
                  (prefix_file_name + "_a_phi.npy"))
                     .string(),
                 a_phi);
    xt::dump_npy((std::filesystem::path{output_dir} / sub_dir /
                  (prefix_file_name + "_f_theta.npy"))
                     .string(),
                 f_theta);
  }
}

auto NFFFTFrequencyDomain::generateSurface() -> void {
  if (!valid()) {
    return;
  }

  for (const auto& f : _frequencies) {
    _fd_plane_data.emplace_back(gridSpace(), emf(), f, nodeTaskSurfaceXN(),
                                nodeTaskSurfaceXP(), nodeTaskSurfaceYN(),
                                nodeTaskSurfaceYP(), nodeTaskSurfaceZN(),
                                nodeTaskSurfaceZP());
  }
}

}  // namespace xfdtd