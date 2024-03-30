#include <xfdtd/network/network.h>
#include <xfdtd/parallel/mpi_support.h>

#include <exception>
#include <filesystem>
#include <utility>
#include <xtensor.hpp>
#include <xtensor/xnpy.hpp>

namespace xfdtd {

Network::Network(std::string output_dir) : _output_dir{std::move(output_dir)} {}

Network::Network(std::vector<std::shared_ptr<Port>> ports,
                 xt::xarray<double> frequencies, std::string output_dir)
    : _ports{std::move(ports)},
      _frequencies{std::move(frequencies)},
      _output_dir{std::move(output_dir)} {}

void Network::addPort(std::shared_ptr<Port> port) {
  _ports.emplace_back(std::move(port));
}

void Network::setFrequencies(xt::xarray<double> frequencies) {
  _frequencies = std::move(frequencies);
}

void Network::setOutputDir(std::string output_dir) {
  _output_dir = std::move(output_dir);
}

void Network::init(
    const std::shared_ptr<const GridSpace>& grid_space,
    const std::shared_ptr<const CalculationParam>& calculation_param,
    const std::shared_ptr<const EMF>& emf) {
  for (auto& port : _ports) {
    port->init(grid_space, calculation_param, emf);
  }

  for (std::size_t i{0}; i < _ports.size(); ++i) {
    if (!_port_set.insert(_ports[i]->index()).second) {
      throw XFDTDNetworkException("Port index must be unique");
    }
  }
}

void Network::output() {
  for (auto& port : _ports) {
    port->calculateSParameters(_frequencies);
  }

  if (!MpiSupport::instance().isRoot()) {
    return;
  }

  for (std::size_t i{0}; i < _ports.size(); ++i) {
    if (!_ports[i]->isSource()) {
      continue;
    }

    for (std::size_t j{0}; j < _ports.size(); ++j) {
      _s_parameters[_ports[j]->index() * 10 + _ports[i]->index()] =
          _ports[j]->b() / _ports[i]->a();
    }
  }
  auto out_path{std::filesystem::path{_output_dir}};

  try {
    if (!std::filesystem::exists(out_path)) {
      std::filesystem::create_directories(out_path);
    }
  } catch (std::exception& e) {
    std::cerr << e.what() << "\n";
    return;
  }

  for (const auto& e : _s_parameters) {
    auto file{out_path / ("s" + std::to_string(e.first) + ".npy")};
    xt::dump_npy(file.string(), e.second);
  }

  auto frequencies_file{out_path / "frequencies.npy"};
  xt::dump_npy(frequencies_file.string(), _frequencies);
}

}  // namespace xfdtd
