#ifndef __XFDTD_CORE_PARALLELIZED_CONFIG_H__
#define __XFDTD_CORE_PARALLELIZED_CONFIG_H__

#include <string>

namespace xfdtd {

class ParallelizedConfigException : public std::exception {
 public:
  static constexpr std::string_view HEADER = "Parallelized Config Exception: ";

  explicit ParallelizedConfigException(const std::string& message)
      : _message{std::string(HEADER) + message} {}

  auto what() const noexcept -> const char* override {
    return _message.c_str();
  }

 private:
  std::string _message;
};

// homogeneous computing environment
class ParallelizedConfig {
 public:
  static constexpr int NUM_DIMS = 3;

  static constexpr int DIM_X_INDEX = 0;

  static constexpr int DIM_Y_INDEX = 1;

  static constexpr int DIM_Z_INDEX = 2;

  ParallelizedConfig() = default;

  ParallelizedConfig(int num_x, int num_y, int num_z, int id = 0, int size = 1,
                     int root = 0);

  auto dims() const -> const int*;

  auto id() const -> int;

  auto size() const -> int;

  auto root() const -> int;

  auto isRoot() const -> bool;

  auto numX() const -> int;

  auto numY() const -> int;

  auto numZ() const -> int;

  virtual auto xPrev() const -> int = 0;

  virtual auto xNext() const -> int = 0;

  virtual auto yPrev() const -> int = 0;

  virtual auto yNext() const -> int = 0;

  virtual auto zPrev() const -> int = 0;

  virtual auto zNext() const -> int = 0;

  virtual auto toString() const -> std::string;

 protected:
  auto doInit() -> void;

  auto setDims(int num_x, int num_y, int num_z) -> void;

  auto setId(int id) -> void;

  auto setSize(int size) -> void;

  auto setRoot(int root) -> void;

 private:
  int _dims[NUM_DIMS]{1, 1, 1};
  int _id{0};
  int _size{1};
  int _root{0};
};

class ThreadConfig : private ParallelizedConfig {
 public:
  static constexpr int ANY_ID = -1;
  static constexpr int INVALID_ID = -2;

  explicit ThreadConfig(int num_x = 1, int num_y = 1, int num_z = 1,
                        int root = 0);

  auto dims() const -> const int*;

  auto size() const -> int;

  auto root() const -> int;

  auto numX() const -> int;

  auto numY() const -> int;

  auto numZ() const -> int;

  auto setId(int id) -> void;

  auto toString() const -> std::string override;

 private:
  auto xPrev() const -> int override;

  auto xNext() const -> int override;

  auto yPrev() const -> int override;

  auto yNext() const -> int override;

  auto zPrev() const -> int override;

  auto zNext() const -> int override;
};

}  // namespace xfdtd

#endif  // __XFDTD_CORE_PARALLELIZED_CONFIG_H__
