#ifndef __XFDTD_CORE_MPI_CONFIG_H__
#define __XFDTD_CORE_MPI_CONFIG_H__

#include <xfdtd/parallel/parallelized_config.h>

#if defined(XFDTD_CORE_WITH_MPI)
#include <mpi.h>
#endif

namespace xfdtd {

class MpiConfig : public ParallelizedConfig {
  friend class MpiSupport;

 public:
  static auto makeSub(const MpiConfig& config, int color, int num_x,
                      int num_y = 1, int num_z = 1, int root = 0) -> MpiConfig;

 public:
  MpiConfig() = default;

  MpiConfig(const MpiConfig&) = delete;

  auto operator=(const MpiConfig&) -> MpiConfig& = delete;

#if defined(XFDTD_CORE_WITH_MPI)
  MpiConfig(MpiConfig&&) noexcept;

  auto operator=(MpiConfig&&) noexcept -> MpiConfig&;
#else
  MpiConfig(MpiConfig&&) noexcept = default;

  auto operator=(MpiConfig&&) noexcept -> MpiConfig& = default;
#endif

  ~MpiConfig();

  auto rank() const -> int;

  auto createCoord() -> void;

  auto xPrev() const -> int override;

  auto xNext() const -> int override;

  auto yPrev() const -> int override;

  auto yNext() const -> int override;

  auto zPrev() const -> int override;

  auto zNext() const -> int override;

  auto toString() const -> std::string override;

#if defined(XFDTD_CORE_WITH_MPI)
  auto comm() const -> MPI_Comm;

  auto cartComm() const -> MPI_Comm;
#endif

 private:
  static auto makeXFDTDComm(int rank, int color, int num_x, int num_y = 1,
                            int num_z = 1, int root = 0) -> MpiConfig;

  int _periods[NUM_DIMS]{0, 0, 0};
  int _reorder{1};

#if defined(XFDTD_CORE_WITH_MPI)
  MPI_Comm _comm{MPI_COMM_NULL};
  MPI_Comm _cart_comm{MPI_COMM_NULL};
  int _x_prev{MPI_PROC_NULL};
  int _x_next{MPI_PROC_NULL};
  int _y_prev{MPI_PROC_NULL};
  int _y_next{MPI_PROC_NULL};
  int _z_prev{MPI_PROC_NULL};
  int _z_next{MPI_PROC_NULL};
#else
  int _x_prev{0};
  int _x_next{0};
  int _y_prev{0};
  int _y_next{0};
  int _z_prev{0};
  int _z_next{0};
#endif
};

}  // namespace xfdtd

#endif  // __XFDTD_CORE_MPI_CONFIG_H__
