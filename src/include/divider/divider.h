#ifndef _XFDTD_LIB_DIVIDER_H_
#define _XFDTD_LIB_DIVIDER_H_

#include <cstddef>
#include <optional>
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
  class Range {
   public:
    Range() = default;
    Range(std::size_t start, std::size_t end) : _start(start), _end(end) {}

    auto operator[](std::size_t index) const {
      if (index == 0) {
        return _start;
      }
      return _end;
    }

    auto& operator[](std::size_t index) {
      if (index == 0) {
        return _start;
      }
      return _end;
    }

    auto operator+(T value) const {
      return Range{_start + value, _end + value};
    }

    auto operator-(T value) const {
      return Range{_start - value, _end - value};
    }

    auto start() const { return _start; }

    auto end() const { return _end; }

    auto size() const { return _end - _start; }

    auto valid() const { return _start < _end; }

    auto toString() const {
      std::stringstream ss;
      ss << "[" << _start << " " << _end << ")";
      return ss.str();
    }

   private:
    std::size_t _start{};
    std::size_t _end{};
  };

  template <typename T>
  struct Task {
    Range<T> _x_range;
    Range<T> _y_range;
    Range<T> _z_range;

    auto toString() const {
      std::stringstream ss;
      ss << "Task: ";
      ss << "x: " << _x_range.toString() << " ";
      ss << "y: " << _y_range.toString() << " ";
      ss << "z: " << _z_range.toString();
      return ss.str();
    }

    auto xRange() const { return _x_range; }

    auto yRange() const { return _y_range; }

    auto zRange() const { return _z_range; }
  };

  using IndexRange = Range<std::size_t>;

  using IndexTask = Task<std::size_t>;

  template <typename T>
  static auto makeRange(T start, T end) {
    return Range<T>{start, end};
  }

  template <typename T>
  static auto makeTask(const Range<T>& x_range, const Range<T>& y_range,
                       const Range<T>& z_range) {
    return Task<T>{x_range, y_range, z_range};
  }

  auto makeIndexRange(std::size_t start, std::size_t end) {
    return makeRange(start, end);
  }

  auto makeIndexTask(const Range<std::size_t>& x_range,
                     const Range<std::size_t>& y_range,
                     const Range<std::size_t>& z_range) {
    return makeTask(x_range, y_range, z_range);
  }

  template <typename T>
  static bool intersected(const Range<T>& domain, const Range<T>& range) {
    return !(range[1] <= domain[0] || domain[1] <= range[0]);
  }

  template <typename T>
  static bool intersected(const Task<T>& task_1, const Task<T>& task_2) {
    return intersected(task_1._x_range, task_2._x_range) &&
           intersected(task_1._y_range, task_2._y_range) &&
           intersected(task_1._z_range, task_2._z_range);
  }

  template <typename T>
  static std::optional<Task<T>> taskIntersection(const Task<T>& task1,
                                                 const Task<T>& task2) {
    if (!intersected(task1, task2)) {
      return {};
    }

    // auto x_range = Range<T>{std::max(task1._x_range[0], task2._x_range[0]),
    //                         std::min(task1._x_range[1], task2._x_range[1])};
    // auto y_range = Range<T>{std::max(task1._y_range[0], task2._y_range[0]),
    //                         std::min(task1._y_range[1], task2._y_range[1])};
    // auto z_range = Range<T>{std::max(task1._z_range[0], task2._z_range[0]),
    //                         std::min(task1._z_range[1], task2._z_range[1])};
    auto x_range =
        Range<T>{std::max(task1.xRange().start(), task2.xRange().start()),
                 std::min(task1.xRange().end(), task2.xRange().end())};
    auto y_range =
        Range<T>{std::max(task1.yRange().start(), task2.yRange().start()),
                 std::min(task1.yRange().end(), task2.yRange().end())};
    auto z_range =
        Range<T>{std::max(task1.zRange().start(), task2.zRange().start()),
                 std::min(task1.zRange().end(), task2.zRange().end())};
    return Task<T>{x_range, y_range, z_range};
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

      res.push_back(makeTask(
          makeRange(x_task.start() + x_start, x_task.end() + x_start),
          makeRange(y_task.start() + y_start, y_task.end() + y_start),
          makeRange(z_task.start() + z_start, z_task.end() + z_start)));
    }

    return res;
  }

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
