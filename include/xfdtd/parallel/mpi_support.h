#ifndef __XFDTD_CORE_MPI_SUPPORT_H__
#define __XFDTD_CORE_MPI_SUPPORT_H__

#include <xfdtd/parallel/mpi_config.h>

#include <complex>
#include <xtensor/xarray.hpp>
#if defined(XFDTD_CORE_WITH_MPI)
#include <mpi.h>
#endif

namespace xfdtd {

class XFDTDMpiSupportException : public XFDTDException {
 public:
  explicit XFDTDMpiSupportException(const std::string& message)
      : XFDTDException("XFDTD Mpi Support Exception: " + message) {}
};

class MpiSupport {
 public:
  class Block {
   public:
    struct Profile {
      int _nx{1};
      int _ny{1};
      int _nz{0};
      int _stride_vec{1};
      int _stride_elem{1};
      int _disp{0};
    };

   public:
    static auto makeBlock(std::size_t nx, std::size_t ny, std::size_t nz)
        -> Block;

    static auto makeExBlock(std::size_t nx, std::size_t ny, std::size_t nz)
        -> Block;

    static auto makeEyBlock(std::size_t nx, std::size_t ny, std::size_t nz)
        -> Block;

    static auto makeEzBlock(std::size_t nx, std::size_t ny, std::size_t nz)
        -> Block;

    static auto makeHxBlock(std::size_t nx, std::size_t ny, std::size_t nz)
        -> Block;

    static auto makeHyBlock(std::size_t nx, std::size_t ny, std::size_t nz)
        -> Block;

    static auto makeHzBlock(std::size_t nx, std::size_t ny, std::size_t nz)
        -> Block;

    static auto make(Profile profile) -> Block;

   public:
    Block() = default;

    Block(Block const&) = delete;

    Block(Block&& other) noexcept;

    auto operator=(Block const&) -> Block& = delete;

    auto operator=(Block&& other) noexcept -> Block&;

    ~Block();

    auto profile() const -> const Profile&;

#if defined(XFDTD_CORE_WITH_MPI)
    auto block() const -> MPI_Datatype { return _block; }
#endif

   private:
    Profile _profile;

#if defined(XFDTD_CORE_WITH_MPI)
    MPI_Datatype _vec_type{MPI_DATATYPE_NULL};

    std::vector<MPI_Datatype> _vec_types;
    std::vector<int> _block_lens;
    std::vector<MPI_Aint> _offsets;

    MPI_Datatype _block{MPI_DATATYPE_NULL};
#endif
  };

  // TODO(franzero): Slice is 2D block. remove it later.
  class Slice {
   public:
    static auto makeRowMajorXSlice(std::size_t nx, std::size_t ny,
                                   std::size_t nz) -> Slice;

    static auto makeRowMajorYSlice(std::size_t nx, std::size_t ny,
                                   std::size_t nz) -> Slice;

    static auto makeRowMajorZSlice(std::size_t nx, std::size_t ny,
                                   std::size_t nz) -> Slice;

    static auto makeEyXSlice(std::size_t nx, std::size_t ny, std::size_t nz)
        -> Slice;

    static auto makeEzXSlice(std::size_t nx, std::size_t ny, std::size_t nz)
        -> Slice;

    static auto makeEzYSlice(std::size_t nx, std::size_t ny, std::size_t nz)
        -> Slice;

    static auto makeExYSlice(std::size_t nx, std::size_t ny, std::size_t nz)
        -> Slice;

    static auto makeExZSlice(std::size_t nx, std::size_t ny, std::size_t nz)
        -> Slice;

    static auto makeEyZSlice(std::size_t nx, std::size_t ny, std::size_t nz)
        -> Slice;

    static auto makeHyXSlice(std::size_t nx, std::size_t ny, std::size_t nz)
        -> Slice;

    static auto makeHzXSlice(std::size_t nx, std::size_t ny, std::size_t nz)
        -> Slice;

    static auto makeHzYSlice(std::size_t nx, std::size_t ny, std::size_t nz)
        -> Slice;

    static auto makeHxYSlice(std::size_t nx, std::size_t ny, std::size_t nz)
        -> Slice;

    static auto makeHxZSlice(std::size_t nx, std::size_t ny, std::size_t nz)
        -> Slice;

    static auto makeHyZSlice(std::size_t nx, std::size_t ny, std::size_t nz)
        -> Slice;

    Slice() = default;

    Slice(Slice const&) = delete;

    Slice(Slice&&) noexcept;

    auto operator=(Slice const&) -> Slice& = delete;

    auto operator=(Slice&&) noexcept -> Slice&;

    ~Slice();

#if defined(XFDTD_CORE_WITH_MPI)
    auto createVec(int count, int block_len, int stripe, MPI_Datatype type)
        -> void;

    auto createStruct(int count, std::vector<int> block_lens,
                      std::vector<MPI_Aint> offsets) -> void;

    auto vecTypes() const -> const std::vector<MPI_Datatype>& {
      return _vec_types;
    }

    auto offsets() const -> const std::vector<MPI_Aint>& { return _offsets; }

    auto blockLens() const -> const std::vector<int>& { return _block_lens; }

    auto slice() const -> MPI_Datatype { return _slice; }
#endif

   private:
#if defined(XFDTD_CORE_WITH_MPI)
    MPI_Datatype _vec_type{MPI_DATATYPE_NULL};

    std::vector<MPI_Datatype> _vec_types;
    std::vector<int> _block_lens;
    std::vector<MPI_Aint> _offsets;

    MPI_Datatype _slice{MPI_DATATYPE_NULL};
#endif
  };

  struct TypeGuard {
    explicit TypeGuard() = default;
    TypeGuard(const TypeGuard&) = delete;
    TypeGuard(TypeGuard&&) = delete;
    TypeGuard& operator=(const TypeGuard&) = delete;
    TypeGuard& operator=(TypeGuard&&) = delete;
    ~TypeGuard() {
#if defined(XFDTD_CORE_WITH_MPI)
      if (_type != MPI_DATATYPE_NULL) {
        MPI_Type_free(&_type);
      }
#endif
    }
#if defined(XFDTD_CORE_WITH_MPI)
    MPI_Datatype _type{MPI_DATATYPE_NULL};
#else
    void* _type{nullptr};
#endif
  };

 public:
#if defined(XFDTD_CORE_WITH_MPI)
  inline static constexpr int ANY_SOURCE = MPI_ANY_SOURCE;
  inline static constexpr int ANY_TAG = MPI_ANY_TAG;
#else
  inline static constexpr int ANY_SOURCE = -1;
  inline static constexpr int ANY_TAG = -1;
#endif

  inline static int exchange_hy_x_sr_tag = 0;
  inline static int exchange_hy_x_rs_tag = 1;
  inline static int exchange_hz_x_sr_tag = 2;
  inline static int exchange_hz_x_rs_tag = 3;

  inline static int exchange_hz_y_sr_tag = 4;
  inline static int exchange_hz_y_rs_tag = 5;
  inline static int exchange_hx_y_sr_tag = 6;
  inline static int exchange_hx_y_rs_tag = 7;

  inline static int exchange_hx_z_sr_tag = 8;
  inline static int exchange_hx_z_rs_tag = 9;
  inline static int exchange_hy_z_sr_tag = 10;
  inline static int exchange_hy_z_rs_tag = 11;

  static auto instance(int argc = 0, char** argv = nullptr) -> MpiSupport& {
    static MpiSupport instance(argc, argv);
    return instance;
  }

 public:
  MpiSupport(MpiSupport const&) = delete;

  MpiSupport(MpiSupport&&) = delete;

  auto operator=(MpiSupport const&) -> MpiSupport& = delete;

  auto operator=(MpiSupport&&) -> MpiSupport& = delete;

  auto serializedPrint(std::string_view message) -> void;

  auto rank() const -> int { return _config.rank(); }

  auto size() const -> int { return _config.size(); }

  auto root() const -> int { return _config.root(); }

  auto isRoot() const -> bool { return _config.root() == _config.rank(); }

  auto abort(int error_code) const -> void;

  auto barrier() const -> void;

  auto generateSlice(std::size_t nx, std::size_t ny, std::size_t nz) -> void;

  auto config() const -> const MpiConfig& { return _config; }

  auto xNext() const -> int;

  auto xPrev() const -> int;

  auto yNext() const -> int;

  auto yPrev() const -> int;

  auto zNext() const -> int;

  auto zPrev() const -> int;

  auto sendRecvHyXHead(xt::xarray<double>& hy) -> void;

  auto recvSendHyXTail(xt::xarray<double>& hy) -> void;

  auto sendRecvHzXHead(xt::xarray<double>& hz) -> void;

  auto recvSendHzXTail(xt::xarray<double>& hz) -> void;

  auto sendRecvHzYHead(xt::xarray<double>& hz) -> void;

  auto recvSendHzYTail(xt::xarray<double>& hz) -> void;

  auto sendRecvHxYHead(xt::xarray<double>& hx) -> void;

  auto recvSendHxYTail(xt::xarray<double>& hx) -> void;

  auto sendRecvHxZHead(xt::xarray<double>& hx) -> void;

  auto recvSendHxZTail(xt::xarray<double>& hx) -> void;

  auto sendRecvHyZHead(xt::xarray<double>& hy) -> void;

  auto recvSendHyZTail(xt::xarray<double>& hy) -> void;

  auto waitAll() -> void;

  auto toString() const -> std::string;

  /**
   * @brief Wrapper for MPI_Send
   *
   * @param config
   * @param buf
   * @param count unit: byte
   * @param dest
   * @param tag
   */
  auto send(const MpiConfig& config, const void* buf, int count, int dest,
            int tag) -> void;

  auto send(const MpiConfig& config, const void* buf, int count,
            const TypeGuard& type, int dest, int tag) -> void;

  auto send(const MpiConfig& config, const double* buf, int count,
            const Block& block, int dest, int tag) -> void;

  /**
   * @brief Wrapper for MPI_Isend. add request to the request list. call waitAll
   * to wait for all requests to be completed. Or call wait(request_id) to wait.
   * Unsafe thread.
   *
   * @param config
   * @param buf
   * @param count
   * @param dest
   * @param tag
   * @return std::size_t: request id
   */
  auto iSend(const MpiConfig& config, const void* buf, int count, int dest,
             int tag) -> std::size_t;

  auto iSend(const MpiConfig& config, const void* buf, int count,
             const TypeGuard& type, int dest, int tag) -> std::size_t;

  auto iSend(const MpiConfig& config, const double* buf, int count,
             const Block& block, int dest, int tag) -> std::size_t;

  auto recv(const MpiConfig& config, void* buf, int count, int source, int tag)
      -> void;

  auto recv(const MpiConfig& config, void* buf, int count,
            const TypeGuard& type, int source, int tag) -> void;

  auto recv(const MpiConfig& config, double* buf, int count, const Block& block,
            int source, int tag) -> void;

  auto iRecv(const MpiConfig& config, void* buf, int count, int source, int tag)
      -> std::size_t;

  auto iRecv(const MpiConfig& config, void* buf, int count,
             const TypeGuard& type, int source, int tag) -> std::size_t;

  auto iRecv(const MpiConfig& config, double* buf, int count,
             const Block& block, int source, int tag) -> std::size_t;

  auto sendRecv(const MpiConfig& config, const void* send_buf, int send_count,
                int dest, int send_tag, void* recv_buf, int recv_count,
                int source, int recv_tag) -> void;

  auto sendRecv(const MpiConfig& config, const void* send_buf, int send_count,
                const TypeGuard& send_type, int dest, int send_tag,
                void* recv_buf, int recv_count, const TypeGuard& recv_type,
                int source, int recv_tag) -> void;

  auto sendRecv(const MpiConfig& config, const double* send_buf, int send_count,
                int dest, int send_tag, double* recv_buf, int recv_count,
                const Block& recv_block, int source, int recv_tag) -> void;

  auto iSendRecv(const MpiConfig& config, const void* send_buf, int send_count,
                 int dest, int send_tag, void* recv_buf, int recv_count,
                 int source, int recv_tag) -> std::size_t;

  auto iSendRecv(const MpiConfig& config, const void* send_buf, int send_count,
                 const TypeGuard& send_type, int dest, int send_tag,
                 void* recv_buf, int recv_count, const TypeGuard& recv_type,
                 int source, int recv_tag) -> std::size_t;

  auto iSendRecv(const MpiConfig& config, const double* send_buf,
                 int send_count, int dest, int send_tag, double* recv_buf,
                 int recv_count, const Block& recv_block, int source,
                 int recv_tag) -> std::size_t;

  auto gather(const MpiConfig& config, const void* send_buf, int send_count,
              void* recv_buf, int recv_count, int root) -> void;

  auto gather(const MpiConfig& config, const void* send_buf, int send_count,
              const TypeGuard& send_type, void* recv_buf, int recv_count,
              const TypeGuard& recv_type, int root) -> void;

  auto gather(const MpiConfig& config, const double* send_buf, int send_count,
              double* recv_buf, int recv_count, const Block& recv_block,
              int root) -> void;

  auto iGather(const MpiConfig& config, const void* send_buf, int send_count,
               void* recv_buf, int recv_count, int root) -> std::size_t;

  auto iGather(const MpiConfig& config, const void* send_buf, int send_count,
               const TypeGuard& send_type, void* recv_buf, int recv_count,
               const TypeGuard& recv_type, int root) -> std::size_t;

  auto iGather(const MpiConfig& config, const double* send_buf, int send_count,
               double* recv_buf, int recv_count, const Block& recv_block,
               int root) -> std::size_t;

  auto allGather(const MpiConfig& config, const void* send_buf, int send_count,
                 void* recv_buf, int recv_count) -> void;

  auto reduceSum(const MpiConfig& config, const double* send_buf,
                 double* recv_buf, int count) -> void;

  auto reduceSum(const MpiConfig& config, const std::complex<double>* send_buf,
                 std::complex<double>* recv_buf, int count) const
      -> void;

 private:
  explicit MpiSupport(int argc = 0, char** argv = nullptr);

  struct MpiGuard {
    ~MpiGuard();
  };
  // this is a guard to finalize MPI when the program ends
  MpiGuard _guard{};

  int _global_rank{0};

  int _global_size{1};

  MpiConfig _config;

  // Slice must be destroyed before MPI_Finalize.
  Slice _hz_y_slice;
  Slice _hx_y_slice;
  Slice _hy_x_slice;
  Slice _hz_x_slice;
  Slice _hx_z_slice;
  Slice _hy_z_slice;
#if defined(XFDTD_CORE_WITH_MPI)
  std::vector<MPI_Request> _requests;
#endif

  auto mpiInit(int argc, char** argv) -> void;
};

}  // namespace xfdtd

#endif  // __XFDTD_CORE_MPI_SUPPORT_H__
