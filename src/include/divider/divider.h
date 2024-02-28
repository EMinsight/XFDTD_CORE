#ifndef _XFDTD_LIB_DIVIDER_H_
#define _XFDTD_LIB_DIVIDER_H_

#include <cstddef>
#include <sstream>
#include <vector>
#include <xtensor/xarray.hpp>
#include <xtensor/xfixed.hpp>

#include "xfdtd/exception/exception.h"

namespace xfdtd {

class XFDTDDividerException : public XFDTDException {
 public:
  explicit XFDTDDividerException(
      const std::string& message = "XFDTD Divider Exception")
      : XFDTDException(message) {}
};

class Divider {
 public:
  enum class Type { UNDEFINED, X, Y, Z, XY, XZ, YZ, XYZ };

  template <typename T>
  using Range = xt::xtensor_fixed<T, xt::xshape<2>>;

  template <typename T>
  struct Task {
    Range<T> _x_range;
    Range<T> _y_range;
    Range<T> _z_range;

    auto toString() {
      std::stringstream ss;
      ss << "Task: ";
      ss << "x: [" << _x_range[0] << " " << _x_range[1] << ") y: ["
         << _y_range[0] << " " << _y_range[1] << ") z: [" << _z_range[0] << " "
         << _z_range[1] << ")";
      return ss.str();
    }
  };

  template <typename T>
  static auto divide(const Task<T>& problem, int num_procs, Type type) {
    switch (type) {
      case Type::X: {
        auto new_problem = problem;
        new_problem._y_range = {0, 1};
        new_problem._z_range = {0, 1};
        auto res = divide(new_problem, num_procs);
        for (auto&& r : res) {
          r._y_range = problem._y_range;
          r._z_range = problem._z_range;
        }
        return res;
      }
      case Type::Y: {
        auto new_problem = problem;
        new_problem._x_range = {0, 1};
        new_problem._z_range = {0, 1};
        auto res = divide(new_problem, num_procs);
        for (auto&& r : res) {
          r._x_range = problem._x_range;
          r._z_range = problem._z_range;
        }
        return res;
      }
      case Type::Z: {
        auto new_problem = problem;
        new_problem._x_range = {0, 1};
        new_problem._y_range = {0, 1};
        auto res = divide(new_problem, num_procs);
        for (auto&& r : res) {
          r._x_range = problem._x_range;
          r._y_range = problem._y_range;
        }
        return res;
      }
      case Type::XY: {
        auto new_problem = problem;
        new_problem._z_range = {0, 1};
        auto res = divide(new_problem, num_procs);
        for (auto&& r : res) {
          r._z_range = problem._z_range;
        }
        return res;
      }
      case Type::XZ: {
        auto new_problem = problem;
        new_problem._y_range = {0, 1};
        auto res = divide(new_problem, num_procs);
        for (auto&& r : res) {
          r._y_range = problem._y_range;
        }
        return res;
      }
      case Type::YZ: {
        auto new_problem = problem;
        new_problem._x_range = {0, 1};
        auto res = divide(new_problem, num_procs);
        for (auto&& r : res) {
          r._x_range = problem._x_range;
        }
        return res;
      }
      case Type::XYZ:
        return divide(problem, num_procs);
      default:
        throw XFDTDDividerException("Invalid type");
    }
  }

  template <typename T>
  static auto divide(const Task<T>& problem, int num_procs) {
    std::vector<Task<std::size_t>> res;

    auto x_start = problem._x_range[0];
    auto x_end = problem._x_range[1];
    auto y_start = problem._y_range[0];
    auto y_end = problem._y_range[1];
    auto z_start = problem._z_range[0];
    auto z_end = problem._z_range[1];

    if (x_end <= x_start || y_end <= y_start || z_end <= z_start) {
      throw XFDTDDividerException("Invalid range");
    }

    auto x_size = x_end - x_start;
    auto y_size = y_end - y_start;
    auto z_size = z_end - z_start;

    for (int i = 0; i < num_procs; ++i) {
      auto x_task = divide(i, num_procs, x_size);
      auto y_task = divide(i, num_procs, y_size);
      auto z_task = divide(i, num_procs, z_size);
      res.push_back({x_task + x_start, y_task + y_start, z_task + z_start});
    }

    return res;
  }

  template <typename T>
  static Range<T> divide(int my_rank, int num_procs, const T& problem_size) {
    if (num_procs <= 0) {
      throw XFDTDDividerException(
          "Number of processes is less than or equal to zero");
    }

    auto quotient = problem_size / num_procs;
    auto remainder = problem_size % num_procs;
    auto start = T{};
    auto end = T{};

    if (my_rank < remainder) {
      start = my_rank * (quotient + 1);
      end = start + quotient + 1;
      return {start, end};
    }

    start = my_rank * quotient + remainder;
    end = start + quotient;
    return {start, end};
  }

  template <typename T>
  auto fromTaskToType(const Task<T>& task) {
    throw XFDTDDividerException("Not implemented");
  }

 private:
};

}  // namespace xfdtd

#endif  // _XFDTD_LIB_DIVIDER_H_
