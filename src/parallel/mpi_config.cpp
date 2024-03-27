#include <xfdtd/parallel/mpi_config.h>

namespace xfdtd {



auto MpiConfig::makeXFDTDComm(int rank, int color, int num_x, int num_y,
                              int num_z, int root) -> MpiConfig {
  MpiConfig mpi_config;
#if defined(XFDTD_CORE_WITH_MPI)
  int size = 1;
  MPI_Comm_split(MPI_COMM_WORLD, color, rank, &mpi_config._comm);
  MPI_Comm_rank(mpi_config._comm, &rank);
  MPI_Comm_size(mpi_config._comm, &size);
  mpi_config.setId(rank);
  mpi_config.setSize(size);
  mpi_config.setDims(num_x, num_y, num_z);
  mpi_config.doInit();
  mpi_config.createCoord();
#endif
  return mpi_config;
}

auto MpiConfig::makeSub(const MpiConfig& config, int color, int num_x,
                        int num_y, int num_z, int root) -> MpiConfig {
  MpiConfig new_config;
#if defined(XFDTD_CORE_WITH_MPI)
  int rank = 0;
  int size = 1;
  MPI_Comm new_comm;
  MPI_Comm_split(config.comm(), color, config.rank(), &new_comm);
  if (new_comm == MPI_COMM_NULL) {
    return new_config;
  }

  new_config._comm = new_comm;

  MPI_Comm_rank(new_config._comm, &rank);
  MPI_Comm_size(new_config._comm, &size);
  new_config.setId(rank);
  new_config.setSize(size);
  new_config.setDims(num_x, num_y, num_z);
  new_config.setRoot(root);
  new_config.doInit();
  new_config.createCoord();
#endif
  return new_config;
}

#if defined(XFDTD_CORE_WITH_MPI)

MpiConfig::MpiConfig(MpiConfig&& other) noexcept
    : ParallelizedConfig(std::move(other)),
      _comm{other._comm},
      _cart_comm{other._cart_comm},
      _x_prev{other._x_prev},
      _x_next{other._x_next},
      _y_prev{other._y_prev},
      _y_next{other._y_next},
      _z_prev{other._z_prev},
      _z_next{other._z_next} {
  other._comm = MPI_COMM_NULL;
  other._cart_comm = MPI_COMM_NULL;
}

auto MpiConfig::operator=(MpiConfig&& other) noexcept -> MpiConfig& {
  if (this != &other) {
    ParallelizedConfig::operator=(std::move(other));
    _comm = other._comm;
    _cart_comm = other._cart_comm;
    _x_prev = other._x_prev;
    _x_next = other._x_next;
    _y_prev = other._y_prev;
    _y_next = other._y_next;
    _z_prev = other._z_prev;
    _z_next = other._z_next;

    other._comm = MPI_COMM_NULL;
    other._cart_comm = MPI_COMM_NULL;
  }

  return *this;
}

MpiConfig::~MpiConfig() {
  if (_cart_comm != MPI_COMM_NULL) {
    MPI_Comm_free(&_cart_comm);
  }

  if (_comm != MPI_COMM_NULL) {
    MPI_Comm_free(&_comm);
  }
}
#else
MpiConfig::~MpiConfig() {}
#endif

auto MpiConfig::rank() const -> int { return id(); }

auto MpiConfig::createCoord() -> void {
#if defined(XFDTD_CORE_WITH_MPI)
  if (_comm == MPI_COMM_NULL) {
    throw ParallelizedConfigException("Invalid MPI communicator");
  }

  MPI_Cart_create(_comm, ParallelizedConfig::NUM_DIMS, dims(), _periods, 1,
                  &_cart_comm);
  MPI_Cart_shift(_cart_comm, ParallelizedConfig::DIM_X_INDEX, 1, &_x_prev,
                 &_x_next);
  MPI_Cart_shift(_cart_comm, ParallelizedConfig::DIM_Y_INDEX, 1, &_y_prev,
                 &_y_next);
  MPI_Cart_shift(_cart_comm, ParallelizedConfig::DIM_Z_INDEX, 1, &_z_prev,
                 &_z_next);
#endif
}

auto MpiConfig::xPrev() const -> int { return _x_prev; }

auto MpiConfig::xNext() const -> int { return _x_next; }

auto MpiConfig::yPrev() const -> int { return _y_prev; }

auto MpiConfig::yNext() const -> int { return _y_next; }

auto MpiConfig::zPrev() const -> int { return _z_prev; }

auto MpiConfig::zNext() const -> int { return _z_next; }

auto MpiConfig::toString() const -> std::string {
  std::stringstream ss;
  ss << "MPI Config:\n";
  ss << " Rank: " << id() << "\n";
  ss << " Num X: " << numX() << "\n";
  ss << " Num Y: " << numY() << "\n";
  ss << " Num Z: " << numZ() << "\n";
  ss << " Size: " << size() << "\n";
  ss << " Root: " << root() << "\n";
  ss << " X: " << xPrev() << "<-" << id() << "->" << xNext() << "\n";
  ss << " Y: " << yPrev() << "<-" << id() << "->" << yNext() << "\n";
  ss << " Z: " << zPrev() << "<-" << id() << "->" << zNext() << "\n";
  ss << " Type: " << Divider::toString(dividerType());
#if defined(XFDTD_CORE_WITH_MPI)
  ss << "\n";
  ss << " MPI Comm: " << comm() << "\n";
  ss << " Cart Comm: " << cartComm();
#endif
  return ss.str();
}

#if defined(XFDTD_CORE_WITH_MPI)
auto MpiConfig::comm() const -> MPI_Comm { return _comm; }

auto MpiConfig::cartComm() const -> MPI_Comm { return _cart_comm; }

#endif

}  // namespace xfdtd