
#include <xfdtd/parallel/mpi_support.h>
#include <xfdtd/parallel/parallelized_config.h>

#include <iostream>
#include <sstream>
#include <string>

#if defined(XFDTD_CORE_WITH_MPI)
#include <mpi.h>
#endif

namespace xfdtd {

ParallelizedConfig::ParallelizedConfig(int num_x, int num_y, int num_z, int id,
                                       int size, int root)
    : _dims{num_x, num_y, num_z}, _id{id}, _size{size}, _root{root} {
  doInit();
}

auto ParallelizedConfig::dims() const -> const int* { return _dims; }

auto ParallelizedConfig::id() const -> int { return _id; }

auto ParallelizedConfig::size() const -> int { return _size; }

auto ParallelizedConfig::root() const -> int { return _root; }

auto ParallelizedConfig::isRoot() const -> bool { return _id == _root; }

auto ParallelizedConfig::numX() const -> int { return _dims[DIM_X_INDEX]; }

auto ParallelizedConfig::numY() const -> int { return _dims[DIM_Y_INDEX]; }

auto ParallelizedConfig::numZ() const -> int { return _dims[DIM_Z_INDEX]; }

auto ParallelizedConfig::toString() const -> std::string {
  std::stringstream ss;
  ss << "Parallelized Config:\n";
  ss << " ID: " << id() << "\n";
  ss << " Num X: " << numX() << "\n";
  ss << " Num Y: " << numY() << "\n";
  ss << " Num Z: " << numZ() << "\n";
  ss << " Size: " << size() << "\n";
  ss << " Root: " << root() << "\n";
  ss << " X: " << xPrev() << "<-" << id() << "->" << xNext() << "\n";
  ss << " Y: " << yPrev() << "<-" << id() << "->" << yNext() << "\n";
  ss << " Z: " << zPrev() << "<-" << id() << "->" << zNext() << "\n";
  return ss.str();
}

auto ParallelizedConfig::doInit() -> void {
  const auto num_x = _dims[DIM_X_INDEX];
  const auto num_y = _dims[DIM_Y_INDEX];
  const auto num_z = _dims[DIM_Z_INDEX];

  if (num_x < 1 || num_y < 1 || num_z < 1) {
    std::stringstream ss;
    ss << "Invalid number of processes: x=" << num_x << ", y=" << num_y
       << ", z=" << num_z;
    throw ParallelizedConfigException(ss.str());
  }

  if (num_x * num_y * num_z != size()) {
    std::stringstream ss;
    ss << "Invalid number of processes: size=" << size();
    throw ParallelizedConfigException(ss.str());
  }

  if (id() < 0 || id() >= size()) {
    std::stringstream ss;
    ss << "Invalid process id: id=" << id();
    throw ParallelizedConfigException(ss.str());
  }

  if (root() < 0 || root() >= size()) {
    std::stringstream ss;
    ss << "Invalid root id: root=" << root();
    throw ParallelizedConfigException(ss.str());
  }
}

auto ParallelizedConfig::setDims(int num_x, int num_y, int num_z) -> void {
  _dims[DIM_X_INDEX] = num_x;
  _dims[DIM_Y_INDEX] = num_y;
  _dims[DIM_Z_INDEX] = num_z;
}

auto ParallelizedConfig::setId(int id) -> void { _id = id; }

auto ParallelizedConfig::setSize(int size) -> void { _size = size; }

auto ParallelizedConfig::setRoot(int root) -> void { _root = root; }

ThreadConfig::ThreadConfig(int num_x, int num_y, int num_z, int root)
    : ParallelizedConfig(num_x, num_y, num_z, 0, num_x * num_y * num_z, root) {}

auto ThreadConfig::dims() const -> const int* {
  return ParallelizedConfig::dims();
}

auto ThreadConfig::size() const -> int { return ParallelizedConfig::size(); }

auto ThreadConfig::root() const -> int { return ParallelizedConfig::root(); }

auto ThreadConfig::numX() const -> int { return ParallelizedConfig::numX(); }

auto ThreadConfig::numY() const -> int { return ParallelizedConfig::numY(); }

auto ThreadConfig::numZ() const -> int { return ParallelizedConfig::numZ(); }

auto ThreadConfig::setId(int id) -> void { ParallelizedConfig::setId(id); }

auto ThreadConfig::toString() const -> std::string {
  std::stringstream ss;
  ss << "ThreadConfig Config:\n";
  ss << " Size: " << size() << "\n";
  ss << " Root: " << root() << "\n";
  return ss.str();
}

auto ThreadConfig::xPrev() const -> int {
  const auto id = this->id();
  const auto ny = numY();
  const auto nz = numZ();

  const auto x = id / (ny * nz);
  return (x == 0) ? (INVALID_ID) : (id - ny * nz);
}

auto ThreadConfig::xNext() const -> int {
  const auto id = this->id();
  const auto ny = numY();
  const auto nz = numZ();

  const auto x = id / (ny * nz);
  return (x == numX() - 1) ? (INVALID_ID) : (id + ny * nz);
}

auto ThreadConfig::yPrev() const -> int {
  const auto id = this->id();
  const auto ny = numY();
  const auto nz = numZ();

  const auto y = (id % (ny * nz)) / nz;
  return (y == 0) ? (INVALID_ID) : (id - nz);
}

auto ThreadConfig::yNext() const -> int {
  const auto id = this->id();
  const auto ny = numY();
  const auto nz = numZ();

  const auto y = (id % (ny * nz)) / nz;
  return (y == ny - 1) ? (INVALID_ID) : (id + nz);
}

auto ThreadConfig::zPrev() const -> int {
  return (id() % numZ() == 0) ? (INVALID_ID) : (id() - 1);
}

auto ThreadConfig::zNext() const -> int {
  return (id() % numZ() == numZ() - 1) ? (INVALID_ID) : (id() + 1);
}

}  // namespace xfdtd
